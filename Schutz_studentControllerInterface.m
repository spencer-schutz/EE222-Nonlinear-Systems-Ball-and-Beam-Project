classdef Schutz_studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        theta_d = 0;

        % S. Schutz - model parameters
        model_params = struct('r_g',0.0254,'L',0.4255,'g',9.81,'K',1.5,'tau',0.025)
        
        % S. Schutz - helper variables
        nx = 4;
        nu = 1;
        
        % S. Schutz - state estimation
        u = 0; % previous input
        x_est = [-.19; 0; 0; 0]; % state estimation
        int_e = 0;

    end

    methods(Access = protected)

        function setupImpl(obj)
            disp("initializing")
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
            %% S. Schutz - General helper variables
            % Extract model parameters
            r_g = obj.model_params.r_g;
            L = obj.model_params.L;
            g = obj.model_params.g;
            K = obj.model_params.K;
            tau = obj.model_params.tau;
            nx = obj.nx;
            nu = obj.nu;

            %% S. Schutz - State Estimator = Cheating w. known X0 + Dynamics
            % Updates property x_est
            t_prev = obj.t_prev;
            x_prev = obj.x_est;
            u_prev = obj.u;
            ode_opt = odeset('Events', @event_ball_out_of_range);
            [t_soln, x_soln] = ode45(@(t,x) ball_and_beam_dynamics(t, x, u_prev), linspace(t_prev, t, 100), x_prev, ode_opt);
            x_est = x_soln(end,:)';
            %disp(x_est)
            obj.x_est = x_est;

            %% S. Schutz - DLQR
            % 1. Linearize SS about x_est
            A_est = [0, 1, 0, 0;
                (5*r_g^2*x_est(4)^2*cos(x_est(3))^2)/(7*L^2), 0, (5*r_g*(2*L*g*cos(x_est(3)) + L*r_g*x_est(4)^2*sin(2*x_est(3)) - 2*r_g*x_est(1)*x_est(4)^2*sin(2*x_est(3))))/(14*L^2), -(5*r_g^2*x_est(4)*cos(x_est(3))^2*(L - 2*x_est(1)))/(7*L^2);
                0, 0, 0, 1;
                0, 0, 0, -1/tau];
            B_est = [0;0;0;K/tau]; 

            % 2. Discretize
            % Euler is ok b/c already assumed const dyn in linearization)
            Ts = t-t_prev;
            A_est_d = eye(nx)+A_est*Ts;
            B_est_d = B_est*Ts;

            % 3. Given r, find x_des and u_des
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            r = p_ball_ref;
            C_r = [1 0 0 0];
            D_r = [0];
            sol = [eye(nx)-A_est_d -B_est_d; C_r D_r]^-1*[zeros(nx,1); eye(nu,1)]*r;
            x_des = sol(1:nx, 1);
            u_des = sol(end,1);

            % 3. Linearize & Discretize (Euler) about x_des
            A_des = [0, 1, 0, 0;
                (5*r_g^2*x_des(4)^2*cos(x_des(3))^2)/(7*L^2), 0, (5*r_g*(2*L*g*cos(x_des(3)) + L*r_g*x_des(4)^2*sin(2*x_des(3)) - 2*r_g*x_des(1)*x_des(4)^2*sin(2*x_des(3))))/(14*L^2), -(5*r_g^2*x_des(4)*cos(x_des(3))^2*(L - 2*x_des(1)))/(7*L^2);
                0, 0, 0, 1;
                0, 0, 0, -1/tau];
            B_des = [0;0;0;K/tau]; 
            A_des_d = eye(nx)+A_des*Ts;
            B_des_d = B_des*Ts;

            % 4. DLQR weights
            Q = diag([2000, 1, 1, 1]);
            R = diag(1);

            % 5. Apply DLQR
            [K_lqr,~,~] = dlqr(A_des_d, B_des_d, Q, R);
            u = u_des - K_lqr*(x_est-x_des);
            V_servo = u;

            % 6. Update properties
            obj.t_prev = t;
            obj.u = u;

            
            %% S. Schutz - Reference Tracking LQR
            % % 1. Linearize SS about x_est
            % A_est = [0, 1, 0, 0;
            %     (5*r_g^2*x_est(4)^2*cos(x_est(3))^2)/(7*L^2), 0, (5*r_g*(2*L*g*cos(x_est(3)) + L*r_g*x_est(4)^2*sin(2*x_est(3)) - 2*r_g*x_est(1)*x_est(4)^2*sin(2*x_est(3))))/(14*L^2), -(5*r_g^2*x_est(4)*cos(x_est(3))^2*(L - 2*x_est(1)))/(7*L^2);
            %     0, 0, 0, 1;
            %     0, 0, 0, -1/tau];
            % B_est = [0;0;0;K/tau]; 
            % 
            % % 2. Given r, find x_des and u_des
            % [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            % r = p_ball_ref;
            % 
            % C_r = [1 0 0 0];
            % D_r = [0];
            % sol = [A_est B_est; C_r D_r]^-1*[zeros(nx,1); eye(nu,1)]*r;
            % x_des = sol(1:nx, 1);
            % u_des = sol(end,1);
            % 
            % % 3. Linearize about x_des
            % A_des = [0, 1, 0, 0;
            %     (5*r_g^2*x_des(4)^2*cos(x_des(3))^2)/(7*L^2), 0, (5*r_g*(2*L*g*cos(x_des(3)) + L*r_g*x_des(4)^2*sin(2*x_des(3)) - 2*r_g*x_des(1)*x_des(4)^2*sin(2*x_des(3))))/(14*L^2), -(5*r_g^2*x_des(4)*cos(x_des(3))^2*(L - 2*x_des(1)))/(7*L^2);
            %     0, 0, 0, 1;
            %     0, 0, 0, -1/tau];
            % B_des = [0;0;0;K/tau]; 
            % 
            % % 4. LQR weights
            % Q = diag([1000, 1, 1, 1]);
            % R = diag(1);
            % 
            % % 5. Apply LQR
            % [K_lqr,~,~] = lqr(A_des, B_des, Q, R);
            % u = u_des - K_lqr*(x_est-x_des);
            % V_servo = u;
            % 
            %  % 6. Update properties
            % obj.t_prev = t;
            % obj.u = u;

            %% S. Schutz - LQI
            % % 1. Linearize SS about x_est
            % A_est = [0, 1, 0, 0;
            %     (5*r_g^2*x_est(4)^2*cos(x_est(3))^2)/(7*L^2), 0, (5*r_g*(2*L*g*cos(x_est(3)) + L*r_g*x_est(4)^2*sin(2*x_est(3)) - 2*r_g*x_est(1)*x_est(4)^2*sin(2*x_est(3))))/(14*L^2), -(5*r_g^2*x_est(4)*cos(x_est(3))^2*(L - 2*x_est(1)))/(7*L^2);
            %     0, 0, 0, 1;
            %     0, 0, 0, -1/tau];
            % B_est = [0;0;0;K/tau]; 
            % 
            % % 2. Get r
            % [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            % r = p_ball_ref;
            % 
            % % 3. Update int_e
            % int_e_old = obj.int_e;
            % dt = (t-t_prev)/99;
            % e = x_soln(:,1)' - r;
            % int_e_since = sum(e*dt);
            % int_e_new = int_e_old + int_e_since;
            % obj.int_e = int_e_new;
            % disp(int_e_new);
            % 
            % % 3. Augmented SS
            % Ai = [0 [1 0 0 0]; zeros(nx, 1) A_est];
            % Bi = [0; B_est];
            % 
            % % 4. LQI Weights
            % Qi = diag([10, zeros(1, nx)]);
            % Ri = 1;
            % 
            % % 5. Apply LQI
            % [K_lqi,~,~] = lqr(Ai, Bi, Qi, Ri);
            % u = -K_lqi*[int_e_new;x_est];
            % V_servo = u;
            % 
            %  % 6. Update properties
            % obj.t_prev = t;
            % obj.u = u; 
            
            %% Sample Controller: Simple Proportional Controller
            % t_prev = obj.t_prev;
            % % Extract reference trajectory at the current timestep.
            % [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            % 
            % % Decide desired servo angle based on simple proportional feedback.
            % k_p = 10; % original value = 3
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
            % k_servo = 5; % original value = 10
            % V_servo = k_servo * (theta_d - theta);
            % obj.u = V_servo;
            % 
            % % Update class properties if necessary.
            % obj.t_prev = t;
            % obj.theta_d = theta_d;
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
