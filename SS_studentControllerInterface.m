classdef studentControllerInterface < matlab.System
    properties (Access = private)
        % General
        t_prev = 0;
        dt_prev = 0.001;
        u = 0;

        % EKF
        xm = [-0.19; 0; 0; 0];
        Pm = diag([0.01, 0.01, 0.01, 0.01]); % P0 = TUNING PARAMETER
    end

    methods(Access = protected)
        %% Step + Step
        function setupImpl(obj)
           disp("Initialize.");
        end

        function [V_servo, xm, dt] = stepImpl(obj, t, p_ball, theta)
            %% 1. Setup
            % coder.extrinsic("dlqr");
            r_g = 0.0254;
            L = 0.4255;
            g = 9.81;
            K = 1.5;
            tau = 0.025;
            
            if t == obj.t_prev
                dt = obj.dt_prev;
            else
                dt = t-obj.t_prev;
            end
            disp([t, dt]) % TODO - delete dt output when done

            [p_ref, v_ref, a_ref] = get_ref_traj(t);

            %% 2. Observer - EKF
            Svv = diag([0.01, 0.01, .01, .01]); % Svv = TUNING PARAMETER
            Sww = diag([0.01, 0.01]); % Sww = TUNING PARAMETER

            % 2a. Predict state - rollout discretized NL dynamics (forward Euler)
            xm_prev = obj.xm;
            u_prev = obj.u;
            dxdt = [xm_prev(2,1);
                    (5*g/7)*(r_g/L)*sin(xm_prev(3,1))-(5/7)*(L/2-xm_prev(1,1))*(r_g/L)^2*xm_prev(4,1)^2*(cos(xm_prev(3,1)))^2;
                    xm_prev(4,1);
                    -xm_prev(4,1)/tau + (K/tau)*u_prev];
            xp_now = xm_prev + dt*dxdt;

            % 2b. Predict variance - linearize discrete NL dyn about xm_prev
            Pm_prev = obj.Pm;
            A_prev =  [1, dt, 0, 0;
                      (5*dt*r_g^2*xm_prev(4,1)^2*cos(xm_prev(3,1))^2)/(7*L^2),  1, (5*dt*r_g*(2*L*g*cos(xm_prev(3,1)) + L*r_g*xm_prev(4,1)^2*sin(2*xm_prev(3,1)) - 2*r_g*xm_prev(1,1)*xm_prev(4,1)^2*sin(2*xm_prev(3,1))))/(14*L^2), -(5*dt*r_g^2*xm_prev(4,1)*cos(xm_prev(3,1))^2*(L - 2*xm_prev(1,1)))/(7*L^2);
                       0,  0,  1, dt;
                       0,  0,  0, 1 - dt/tau]; 
            L_prev = eye(4);
            Pp_now = A_prev*Pm_prev*A_prev' + L_prev*Svv*L_prev';

            % 2c. Measurement update - estimate and variance
            H_now = [1 0 0 0; 0 0 1 0];
            M_now = eye(2);
            z_now = [p_ball; theta];
            zp_now = H_now*xp_now;
            K_now = Pp_now*H_now'*(H_now*Pp_now*H_now' + M_now*Sww*M_now')^-1;
            xm_now = xp_now + K_now*(z_now - zp_now);
            Pm_now = (eye(4) - K_now*H_now) * Pp_now * (eye(4)-K_now*H_now)' + K_now*M_now*Sww*M_now'*K_now';

            % 2d. Update Parameters
            obj.xm = xm_now;
            obj.Pm = Pm_now;
           
            %% 3. Local Linearization-Based LQR
            % 3a. Linearize discrete NL dynamics about desired (x_star, u_star)
            x_star = [p_ref; v_ref; 0; 0];
            u_star = 0;
            A_star =  [1, dt, 0, 0;
                       (5*dt*r_g^2*x_star(4,1)^2*cos(x_star(3,1))^2)/(7*L^2),  1, (5*dt*r_g*(2*L*g*cos(x_star(3,1)) + L*r_g*x_star(4,1)^2*sin(2*x_star(3,1)) - 2*r_g*x_star(1,1)*x_star(4,1)^2*sin(2*x_star(3,1))))/(14*L^2), -(5*dt*r_g^2*x_star(4,1)*cos(x_star(3,1))^2*(L - 2*x_star(1,1)))/(7*L^2);
                       0,  0,  1, dt;
                       0,  0,  0, 1 - dt/tau]; 
            B_star = [0;0;0;dt*(K/tau)];

            % 3b. LQR Weights
            Q = diag([1000, 500, 1, 1]);
            R = diag(0.1);
            
            % 3c. DARE
            max_iter = 100000;
            P = eye(4);
            K = zeros(1,4);
            for i = 1:max_iter
                P_new = A_star'*P*A_star - (A_star'*P*B_star)*(R+B_star'*P*B_star)^-1*(B_star'*P*A_star) + Q;
                K_new = (R+B_star'*P*B_star)^-1*B_star'*P*A_star;
                if isequal(K, K_new)
                    %j = int64(i);
                    %fprintf("Converged on Iteration %d\n", j);
                    break
                elseif i == max_iter
                    j = int64(i);
                    fprintf("Failed to Converge: max_iter = %d\n", j);
                end
                P = P_new;
                K = K_new;
            end
            
            % 3c. Apply LQR
            u_now = u_star - K*(xm_now-x_star);

            % 3d. Update Parameters
            obj.u = u_now; 

            % disp(A_star)
            % disp(B_star)

            % Check Ctrb, Obsv
            % disp(rank([B_star A_star*B_star A_star^2*B_star A_star^3*B_star]));
            % C = [1 0 0 0; 0 1 0 0];
            % disp(rank([C; C*A_star; C*A_star^2; C*A_star^3]));

            %% End: Function Outputs, Step to Next Time
            V_servo = u_now;
            xm = xm_now; % TODO - DELETE EXTRA OUTPUT WHEN DONE 
            obj.t_prev = t;

        end

    end
    
    methods(Access = public)
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)        
            % Run setup at start
            if t==0
                setupImpl(obj);
            end
            % Step controller, observer, etc.
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
    end
    
end