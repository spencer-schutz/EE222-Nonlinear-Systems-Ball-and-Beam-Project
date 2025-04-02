classdef SY_studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and update while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        theta_d = 0;
        
        % S. Schutz - model parameters
        model_params = struct('r_g',0.0254,'L',0.4255,'g',9.81,'K',1.5,'tau',0.025)

        nx = 4;
        nu = 1;
        
        u = 0;
        x_est = [-0.19; 0; 0; 0];
        sig_x = 0; sig_theta = 0;

        d2traj_ref = [];

        
    end
    methods(Access = protected)
        function setupImpl(obj)
           disp("You can use this function for initializaition.");
        end

        function V_servo = stepImpl(obj, t, p_ball, theta)
        % This is the main function called every iteration. You have to implement
        % the controller in this function, bu you are not allowed to
        % change the signature of this function. 
        % Input arguments:
        %   t: current time
        %   p_ball: position of the ball provided by the ball position sensor (m)
        %
        %   theta: servo motor angle provided by the encoder of the motor (rad)
        % Output:
        %   V_servo: voltage to the servo input.        
            %% Debugging
            %tic;
        
            %% Parameters
            r_g = obj.model_params.r_g;
            L = obj.model_params.L;
            g = obj.model_params.g;
            K = obj.model_params.K;
            tau = obj.model_params.tau;
            nx = obj.nx;
            nu = obj.nu;

            sat_sig_x = 0.8;
            sat_sig_theta = 0.8;

            a = 5/7*g*r_g/L;

            %% Update observer (Extended High Gain Observer)
            % Gains
            epsilon = 0.04;
            
            % State parameters
            t_prev = obj.t_prev; 
            dt = t - t_prev;
            x_prev = obj.x_est;
            sig_x_prev = obj.sig_x; sig_theta_prev = obj.sig_theta;
            u_prev = obj.u;
            
            % Observer update
            x_est = zeros(4,1);
            x_est(1) = x_prev(1) + dt * (x_prev(2) + 3/epsilon*(p_ball - x_prev(1)));
            x_est(2) = x_prev(2) + dt * (a * sin(x_prev(1)) + sig_theta_prev + 3/epsilon^2*(p_ball - x_prev(1)));
            sig_x = sig_x_prev + dt * 1/epsilon^3*(p_ball - x_prev(1));
            
            x_est(3) = x_prev(3) + dt * (x_prev(4) + 3/epsilon*(theta - x_prev(3)));
            x_est(4) = x_prev(4) + dt * (-x_prev(4)/tau + K/tau * u_prev + sig_theta_prev + 9/epsilon^2*(theta - x_prev(3)));
            sig_theta = sig_theta_prev + dt * 27/epsilon^3*(theta - x_prev(3));
            
            sig_x = min(abs(sig_x),sat_sig_x)*sign(sig_x);
            sig_theta = min(abs(sig_theta),sat_sig_theta)*sign(sig_theta);
            
            % Update properites
            obj.x_est = x_est;
            obj.sig_x = sig_x; obj.sig_theta = sig_theta;
            
            %% Fit polynomial to get derivative
            [p_ref, dp_ref, d2p_ref] = get_ref_traj(t);
            d2traj_ref = [obj.d2traj_ref, d2p_ref];
            windowSize = 7;
            polyOrder = 2;

            if t >= 0.1
                p = polyfit(t-dt*windowSize:dt:t, d2traj_ref(end-windowSize:end), polyOrder);
                d3p_ref = polyval(polyder(p), t);
                d4p_ref = polyval(polyder(polyder(p)), t);
            else
                d3p_ref = 0;
                d4p_ref = 0;
            end
            obj.d2traj_ref = d2traj_ref;

            %% Controller: LQR + Approximate Feedback Linearization
            % negative pole st. charac. eq. of A = (s - lambda)^4
            % lambda = -3; 
            % k1 = lambda^4;
            % k2 = -4*lambda^3;
            % k3 = 6*lambda^2;
            % k4 = -4*lambda;

            % Get reference traj and calculate third, fourth order derivative
            [p_ref, dp_ref, d2p_ref] = get_ref_traj(t);
            %CHEATING
            % amplitude = 0.04;
            % period = 10; % sec
            % omega = 2 * pi / period; 
            % %a_ref = - amplitude * omega^2 * sin(omega * t);
            % d3p_ref = - amplitude * omega^3 * cos(omega * t);
            % d4p_ref = amplitude * omega^4 * sin(omega * t);

            % Approximate Feedback Linearization and DLQR setup
            f_x = -a/tau * x_est(4) * cos(x_est(3)) - a * x_est(4)^2 * sin(x_est(3));
            g_x = a*K/tau * cos(x_est(3));

            A = [0 1 0 0;
                 0 0 1 0;
                 0 0 0 1;
                 0 0 0 0];
            B = [0 0 0 1]';
            % discretize - might want to pull up to initialization (does dt change?)
            sys = ss(A,B,[1 0 0 0], 0);
            sysd = c2d(sys, dt);
            Ad = sysd.A;
            Bd = sysd.B;
            
            % Apply DLQR
            Q = diag([800, 0.01, 1, 1]);
            R = diag(0.5);
            [K_lqr, ~, ~] = dlqr(Ad, Bd, Q, R);

            % Get control
            e1 = x_est(1) - p_ref;
            e2 = x_est(2) - dp_ref;
            e3 = a*sin(x_est(3)) - d2p_ref;
            e4 = a*x_est(4)*cos(x_est(3)) - d3p_ref;
            V_servo = 1/g_x * (-f_x + d4p_ref - K_lqr(1)*e1 - K_lqr(2)*e2 - K_lqr(3)*e3 - K_lqr(4)*e4);

            % Update class properites if necessary
            obj.t_prev = t;
            obj.u = V_servo;
            
            %elapsed_time = toc;
            %fprintf("elapsed time: %f\n", elapsed_time);
            %% Sample Controller: Simple Proportional Controller
            % t_prev = obj.t_prev;
            % % Extract reference trajectory at the current timestep.
            % [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            % % Decide desired servo angle based on simple proportional feedback.
            % k_p = 3;
            % theta_d = - k_p * (p_ball - p_ball_ref);
            % 
            % % Make sure that the desired servo angle does not exceed the physical
            % % limit. This part of code is not necessary but highly recommended
            % % because it addresses the actual physical limit of the servo motor.
            % theta_saturation = 56 * pi / 180;    
            % theta_d = min(theta_d, theta_saturation);
            % theta_d = max(theta_d, -theta_saturation);
            % 
            % % Simple position control to control servo angle to the desired
            % % position.
            % k_servo = 10;
            % V_servo = k_servo * (theta_d - theta);
            % 
            % % Update class properties if necessary.
            % obj.t_prev = t;
            % obj.theta_d = theta_d;
        end
    end
    
    methods(Access = public)
        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)   
            if t == 0
                setupImpl(obj);
            end
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
    end
    
end