classdef Bosio_studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        theta_d = 0;
        
        r_g = 0.0254;
        L = 0.4255;
        g = 9.81;
        K = 1.5;
        tau = 0.025

        nx = 4;
        nu = 1;

        u = 0; % previous input

        x_est = [0; 0; 0; 0]; % state estimate, and initial condition
        int_e = 0;

        P_m = ones(4,4); % initial state covariance

        sigma_v = 0.1*eye(4); % process noise covariance
        sigma_w = 0.1*eye(2); % measurement noise covariance

        theta_saturation = 56 * pi / 180;
        
        % dynamics function
        f = @(x, u, r_g, L, g, K, tau) ...
                 [x(2);
                  (5*g/7)*(r_g/L)*sin(x(3))-(5/7)*(L/2-x(1))*(r_g/L)^2*x(4)^2*(cos(x(3)))^2;
                  x(4);
                  -x(4)/tau + (K/tau)*u];

        A = @(x, r_g, L, g, K, tau) ...
               [0, 1, 0, 0;
                (5*r_g^2*x(4)^2*cos(x(3))^2)/(7*L^2), 0, (5*r_g*(2*L*g*cos(x(3)) + L*r_g*x(4)^2*sin(2*x(3)) - 2*r_g*x(1)*x(4)^2*sin(2*x(3))))/(14*L^2), -(5*r_g^2*x(4)*cos(x(3))^2*(L - 2*x(1)))/(7*L^2);
                0, 0, 0, 1;
                0, 0, 0, -1/tau];

        H = [1, 0, 0, 0; % measurement matrix
             0, 0, 1, 0];
    end
    methods(Access = protected)
        function setupImpl(obj)
           disp("Initialization.");
        end

        function V_servo = stepImpl(obj, t, p_ball, theta)
        % This is the main function called every iteration. You have to implement
        % the controller in this function, but you are not allowed to
        % change the signature of this function. 
        % Input arguments:
        %   t: current time
        %   p_ball: position of the ball provided by the ball position sensor (m)
        %
        %   theta: servo motor angle provided by the encoder of the motor (rad)
        % Output:
        %   V_servo: voltage to the servo input.        
            %% Sample Controller: Simple Proportional Controller
            t_prev = obj.t_prev;
            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);

            p_est = Estimator(obj, t, p_ball, theta);
            % Decide desired servo angle based on simple proportional feedback.
            k_p = 3;
            k_d = 10;

            obj.theta_d = - k_p * (p_est - p_ball_ref) - k_d*obj.x_est(2);

            % Make sure that the desired servo angle does not exceed the physical
            % limit. This part of code is not necessary but highly recommended
            % because it addresses the actual physical limit of the servo motor.
                
            obj.theta_d = min(obj.theta_d, obj.theta_saturation);
            obj.theta_d = max(obj.theta_d, -obj.theta_saturation);

            % Simple position control to control servo angle to the desired
            % position.
            k_servo = 5;
            V_servo = k_servo * (obj.theta_d - theta);
            
            % Update class properties if necessary.
            obj.t_prev = t;
        end

        function p_est = Estimator(obj, t, p_ball, theta)

            dt = t - obj.t_prev;
             
            % prediction step
            x_p = obj.f(obj.x_est, obj.u, obj.r_g, obj.L, obj.g, obj.K, obj.tau);
            
            % linearization
            A = obj.A(obj.x_est, obj.r_g, obj.L, obj.g, obj.K, obj.tau);
            X = eye(4,4);
            P_p = A*obj.P_m*A' + X*obj.sigma_v*X';

            % measurement update
            Kk = P_p*obj.H'*inv(obj.H*P_p*obj.H' + obj.sigma_w);
            z = [p_ball; theta];
            obj.x_est = x_p + Kk*(z - obj.H*x_p);
            obj.P_m = (eye(4) - Kk*obj.H)*P_p;
            
            p_est = obj.x_est(1);
        end
    end
    
    methods(Access = public)
        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function initController(obj)
            setupImpl(obj)
        end

        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)        
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
    end
    
end