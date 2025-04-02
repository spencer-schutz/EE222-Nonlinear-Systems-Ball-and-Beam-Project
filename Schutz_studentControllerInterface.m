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
        nx = 4; % number of states
        nz = 2; % number of measurements
        nu = 1; % number of inputs
        
        % S. Schutz - state estimation
        u = 0; % init with 0
        xm = [-.19; 0; 0; 0]; % init with x0
        Pm = diag([0.05, 0, 0.01, 0]); % init with P0

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
            nz = obj.nz;
            nu = obj.nu;
            
            %% S. Schutz - State Estimation = EKF (DT)
            % 0. Design Parameters - process and measurement variances
            Svv = diag([0.05, 0.001, .0175, .00175]);
            Sww = diag([0.05, 0.0175]);

            % 1. Predict state - rollout discretized NL dynamics (forward Euler)
            xm_prev = obj.xm; % from previous iteration 
            u_prev = obj.u; % known
            dxdt = [xm_prev(2);
                    (5*g/7)*(r_g/L)*sin(xm_prev(3)) - (5/7) * (L/2-xm_prev(1)) * (r_g/L)^2 * xm_prev(4)^2 * (cos(xm_prev(3)))^2 ;
                    xm_prev(4);
                    -xm_prev(4)/tau + (K/tau)*u_prev]; % CT NL dynamics
            dt = t - obj.t_prev; % get sample time
            xp_now = xm_prev + dt*dxdt; % forward Euler

            % 2. Predict variance - linearize about last estimate
            Pm_prev = obj.Pm; % from previous iteration 
            A_prev = [1, dt, 0, 0;
                    (5*dt*r_g^2*xm_prev(4)^2*cos(xm_prev(3))^2)/(7*L^2),  1, (5*dt*r_g*(2*L*g*cos(xm_prev(3)) + L*r_g*xm_prev(4)^2*sin(2*xm_prev(3)) - 2*r_g*xm_prev(1)*xm_prev(4)^2*sin(2*xm_prev(3))))/(14*L^2), -(5*dt*r_g^2*xm_prev(4)*cos(xm_prev(3))^2*(L - 2*xm_prev(1)))/(7*L^2);
                     0,  0,  1, dt;
                     0,  0,  0, 1 - dt/tau]; % linearize DT NL dyn about xm_prev
            L_prev = dt*eye(nx); % linearize DT NL dyn about v
            Pp_now = A_prev*Pm_prev*A_prev' + L_prev*Svv*L_prev';

            % 3. Measurement update
            H_now = [1 0 0 0; 0 0 1 0];
            M_now = eye(nz);
            z_now = [p_ball; theta];
            zp_now = H_now*xp_now;
            K_now = Pp_now*H_now'*(H_now*Pp_now*H_now' + M_now*Sww*M_now')^-1;
            xm_now = xp_now + K_now*(z_now - zp_now);
            Pm_now = (eye(nx) - K_now*H_now) * Pp_now * (eye(nx)-K_now*H_now)' + K_now*M_now*Sww*M_now'*K_now';

            % Save new values
            obj.xm = xm_now;
            obj.Pm = Pm_now;

            %% S. Schutz - Controller = Reference Tracking LQR (DT)
            % 1. Get reference
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            x_star = [p_ball_ref; v_ball_ref; 0; 0];
            u_star = 0;

            % 2. Linearize DT Dyn about current reference
            A_star = [1, dt, 0, 0;
                    (5*dt*r_g^2*x_star(4)^2*cos(x_star(3))^2)/(7*L^2),  1, (5*dt*r_g*(2*L*g*cos(x_star(3)) + L*r_g*x_star(4)^2*sin(2*x_star(3)) - 2*r_g*x_star(1)*x_star(4)^2*sin(2*x_star(3))))/(14*L^2), -(5*dt*r_g^2*x_star(4)*cos(x_star(3))^2*(L - 2*x_star(1)))/(7*L^2);
                     0,  0,  1, dt;
                     0,  0,  0, 1 - dt/tau]; % linearize DT NL dyn about xm_prev
            B_star = [0;0;0;dt*(K/tau)];

            % 4. LQR weights
            Q = diag([1000, 500, 1, 1]);
            R = diag(0.1);

            % 5. Apply LQR
            [K_lqr,~,~] = dlqr(A_star, B_star, Q, R);
            u = u_star - K_lqr*(xm_now-x_star);
            V_servo = u;

            % 6. Update properties
            obj.t_prev = t;
            obj.u = u;

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
