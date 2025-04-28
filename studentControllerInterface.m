classdef studentControllerInterface < matlab.System
    properties (Access = private)
        % General
        t_prev = 0;
        dt_prev = 0.001;
        u = 0;

        % EKF
        xm = [0; 0; 0; 0];
        Pm = diag([0.01, 0.01, 0.01, 0.01]); % P0 = TUNING PARAMETER

        % Integral
%         sum_e = -0.4; % for LQR
        sum_e = 0; % for FL

    end

    methods(Access = protected)
        %% Step + Step
        function setupImpl(obj)
           disp("Initialize.");
        end

        function [V_servo, xm, xm_dot, theta, theta_dot, sum_e] = stepImpl(obj, t, p_ball, theta)
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
%             disp([t, dt]) % TODO - delete dt output when done

            [p_ref, v_ref, a_ref] = get_ref_traj(t);

            %% 2. Observer - EKF
            Svv = diag([0.01, 0.05, .02, .2]);
            Sww = diag([1, 0.1]);

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
%             % 3a. Linearize discrete NL dynamics about desired (x_star, u_star)
%             x_star = [p_ref; v_ref; 0; 0];
%             u_star = 0;
%             A_star =  [1, dt, 0, 0;
%                        (5*dt*r_g^2*x_star(4,1)^2*cos(x_star(3,1))^2)/(7*L^2),  1, (5*dt*r_g*(2*L*g*cos(x_star(3,1)) + L*r_g*x_star(4,1)^2*sin(2*x_star(3,1)) - 2*r_g*x_star(1,1)*x_star(4,1)^2*sin(2*x_star(3,1))))/(14*L^2), -(5*dt*r_g^2*x_star(4,1)*cos(x_star(3,1))^2*(L - 2*x_star(1,1)))/(7*L^2);
%                        0,  0,  1, dt;
%                        0,  0,  0, 1 - dt/tau]; 
%             B_star = [0;0;0;dt*(K/tau)];
% 
%             % 3b. LQR Weights
%             Q = diag([1000, 100, 0.1, 1e-10]);
%             R = diag(0.4);
%             
%             % 3c. DARE
%             max_iter = 100000;
%             P = eye(4);
%             K_lqr = zeros(1,4);
%             for i = 1:max_iter
%                 P_new = A_star'*P*A_star - (A_star'*P*B_star)*(R+B_star'*P*B_star)^-1*(B_star'*P*A_star) + Q;
%                 K_new = (R+B_star'*P*B_star)^-1*B_star'*P*A_star;
%                 if isequal(K_lqr, K_new)
%                     %j = int64(i);
%                     %fprintf("Converged on Iteration %d\n", j);
%                     break
%                 elseif i == max_iter
%                     j = int32(i);
%                     fprintf("Failed to Converge: max_iter = %d\n", j);
%                 end
%                 P = P_new;
%                 K_lqr = K_new;
%             end
%             
%             % 3c. Apply LQR
%             u_now = u_star - K_lqr*(xm_now-x_star);
%             
%             % 3.c.1. Add integral
% %             Ke = 3;
%             Ke = 0;
%             e_now = xm_now(1) - p_ref;
%             obj.sum_e = obj.sum_e + e_now*dt;
%             obj.sum_e = max(min(obj.sum_e, 1), -1);
%             u_e = Ke*obj.sum_e;
%             u_now = u_now - u_e;
% 
%             % 3d. Update Parameters
%             obj.u = u_now; 
            
            %% FL+LQR
            a = 5/7*g*r_g/L;

            % Approx feedback linearization
            f_x = -a/tau * xm_now(4) * cos(xm_now(3)) - a * xm_now(4)^2 * sin(xm_now(3));
            g_x = a*K/tau * cos(xm_now(3));

            A = [0 1 0 0;
                 0 0 1 0;
                 0 0 0 1;
                 0 0 0 0];
            B = [0 0 0 1]';
            % discretize - might want to pull up to initialization (does dt change?)
            Ad = expm(A*dt);
            n = size(A,1);
            
            % Check if A is invertible
            if rank(A) == n
                % Use closed-form formula for Bd
                Bd = A \ (Ad - eye(n)) * B;
            else
                % A is singular; perform numerical integration
                Nsteps = 100;  % Number of integration steps
                tau = linspace(0, dt, Nsteps);
                dtau = dt/(Nsteps - 1);
                Bd = zeros(n, size(B,2));
                for i = 1:Nsteps
                    Bd = Bd + expm(A * tau(i)) * B * dtau;
                end
            end

            % 3b. LQR Weights
            Q = diag([1000, 100, 0.1, 1e-10]);
            R = diag(0.004);
            
            % 3c. DARE
            max_iter = 100000;
            P = eye(4);
            K_lqr = zeros(1,4);
            for i = 1:max_iter
                P_new = Ad'*P*Ad - (Ad'*P*Bd)*(R+Bd'*P*Bd)^-1*(Bd'*P*Ad) + Q;
                K_new = (R+Bd'*P*Bd)^-1*Bd'*P*Ad;
                if isequal(K_lqr, K_new)
                    %j = int64(i);
                    %fprintf("Converged on Iteration %d\n", j);
                    break
                elseif i == max_iter
                    j = int32(i);
                    fprintf("Failed to Converge: max_iter = %d\n", j);
                end
                P = P_new;
                K_lqr = K_new;
            end

            e1 = xm_now(1) - p_ref;
            e2 = xm_now(2) - v_ref;
            e3 = a*sin(xm_now(3)) - a_ref;
            e4 = a*xm_now(4)*cos(xm_now(3)) - 0;
            u_now = 1/g_x * (-f_x + 0 - K_lqr(1)*e1 - K_lqr(2)*e2 - K_lqr(3)*e3 - K_lqr(4)*e4);

            % 3.c.1. Add integral
            Ke = 2;
            e_now = xm_now(1) - p_ref;
            obj.sum_e = obj.sum_e + e_now*dt;
            obj.sum_e = max(min(obj.sum_e, 1), -1);
            u_e = Ke*obj.sum_e;
            u_now = u_now - u_e;

            % Update parameters
            obj.u = u_now;
            
            %% End: Function Outputs, Step to Next Time
            V_servo = u_now;
            xm = xm_now(1,1); % TODO - DELETE EXTRA OUTPUT WHEN DONE 
            xm_dot = xm_now(2,1);
            theta = xm_now(3,1);
            theta_dot = xm_now(4,1);
            obj.t_prev = t;
            sum_e = obj.sum_e;

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