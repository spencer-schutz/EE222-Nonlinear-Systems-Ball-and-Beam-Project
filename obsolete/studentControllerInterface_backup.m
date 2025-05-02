classdef studentControllerInterface < matlab.System
    properties (Access = private)
        % Time, reference, measurement, input
        t = 0;
        t_prev = -1;
        dt = 0;
        ref = 0;
        p_ball = 0;
        theta = 0;
        V_servo = 0;   

        % Measurement Equation (for EKF)
        H_d = [1 0 0 0; 0 0 1 0];
        M_d = eye(2);

        % EKF variables
        xm = [-0.19; 0; 0; 0]; % init with x0
        Pm = diag([0.05, 0, 0.01, 0]); % init with P0

        % Feedback Linearization variables
        sig_x = 0;
        sig_theta = 0;
        %d2traj_ref = 0;

        theta_d = 0;        
        theta_saturation = 56 * pi / 180;

    end

    methods(Access = protected)
        %% Step + Step
        function setupImpl(obj)
           disp("You can use this function for initializaition.");
        end

        function V_servo = stepImpl(obj, t, p_ball, theta)
            % Publish current step info
            [p_ref, v_ref, a_ref] = get_ref_traj(t);
            
            obj.p_ball = p_ball;
            obj.theta = theta;
            
            obj.t = t;
            obj.dt = t - obj.t_prev;

            % Observer
            EKF(obj)
            
            % Choose controller 
            % Option 1: Switching
            if abs(obj.xm(1) - p_ref) >= 0.02
                LQR(obj)
            else
                FL_LQR(obj)
            end

            % Option 2: approx feedback linearization + LQR
            % FL_LQR(obj)
            
            % Option 3: LQR
            % LQR(obj)
            
            % Save output
            V_servo = obj.V_servo;

            % Go to next step
            obj.t_prev = t;
        end

        %% Observers
        function EKF(obj)
            % 0. Setup
            r_g = 0.0254;
            L = 0.4255;
            g = 9.81;
            K = 1.5;
            tau = 0.025;

            x_est = obj.xm;
            V_servo = obj.V_servo;
            dt = obj.dt;

            Svv = diag([0.05, 0.001, .0175, .00175]);
            Sww = diag([0.05, 0.0175]);
             
            % prediction step
            f_c = [x_est(2);
                  (5*g/7)*(r_g/L)*sin(x_est(3))-(5/7)*(L/2-x_est(1))*(r_g/L)^2*x_est(4)^2*(cos(x_est(3)))^2;
                  x_est(4);
                  -x_est(4)/tau + (K/tau)*V_servo];
            xp = x_est + dt*f_c;
            
            % linearization
            A = [0, 1, 0, 0;
                (5*r_g^2*x_est(4)^2*cos(x_est(3))^2)/(7*L^2), 0, (5*r_g*(2*L*g*cos(x_est(3)) + L*r_g*x_est(4)^2*sin(2*x_est(3)) - 2*r_g*x_est(1)*x_est(4)^2*sin(2*x_est(3))))/(14*L^2), -(5*r_g^2*x_est(4)*cos(x_est(3))^2*(L - 2*x_est(1)))/(7*L^2);
                0, 0, 0, 1;
                0, 0, 0, -1/tau];
            
            X = eye(4,4);
            Pp = A*obj.Pm*A' + X*Svv*X';

            % measurement update
            H = obj.H_d;
            Kk = Pp*H'*inv(H*Pp*H' + Sww);
            z = [obj.p_ball; obj.theta];
            obj.xm = xp + Kk*(z - H*xp);
            obj.Pm = (eye(4) - Kk*H)*Pp;
            
        end

        %% Controllers
        function FL_LQR(obj)
            % 0. Setup
            r_g = 0.0254;
            L = 0.4255;
            g = 9.81;
            K = 1.5;
            tau = 0.025;

            a = 5/7*g*r_g/L;

            % Fit polynomial to get derivative
            [p_ref, dp_ref, d2p_ref] = get_ref_traj(obj.t);
            % ref = obj.ref;
            % p_ref = ref(1);
            % dp_ref = ref(2);
            % d2p_ref = ref(3);
            %d2traj_ref = [obj.d2traj_ref, d2p_ref];
            windowSize = 7;
            polyOrder = 2;
            
            t = obj.t;
            dt = obj.dt;
            
            % Get higher order derivatives of reference
            if t >= 0.1
                %p = polyfit(t-dt*windowSize:dt:t, d2traj_ref(end-windowSize:end), polyOrder);
                d3p_ref = 0; %polyval(polyder(p), t);
                d4p_ref = 0; %polyval(polyder(polyder(p)), t);
            else
                d3p_ref = 0;
                d4p_ref = 0;
            end
            %obj.d2traj_ref = d2traj_ref;

            % Approximate Feedback Linearization and DLQR setup
            x_est = obj.xm;
            f_x = -a/tau * x_est(4) * cos(x_est(3)) - a * x_est(4)^2 * sin(x_est(3));
            g_x = a*K/tau * cos(x_est(3));

            A = [0 1 0 0;
                 0 0 1 0;
                 0 0 0 1;
                 0 0 0 0];
            B = [0 0 0 1]';
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
            
            % Apply DLQR
            Q = diag([800, 50, 1, 1]);
            R = diag(0.3);
            % [K_lqr, ~, ~] = dlqr(Ad, Bd, Q, R);
            K_lqr = dlqr_iterative(obj, Ad, Bd, Q, R);

            % Get control
            e1 = x_est(1) - p_ref;
            e2 = x_est(2) - dp_ref;
            e3 = a*sin(x_est(3)) - d2p_ref;
            e4 = a*x_est(4)*cos(x_est(3)) - d3p_ref;
            V_servo = 1/g_x * (-f_x + obj.sig_theta + d4p_ref - K_lqr(1)*e1 - K_lqr(2)*e2 - K_lqr(3)*e3 - K_lqr(4)*e4);

            % Update properties
            obj.V_servo = V_servo;
        end

        function LQR(obj)
            % 0. Constants
            r_g = 0.0254;
            L = 0.4255;
            g = 9.81;
            K = 1.5;
            tau = 0.025;

            dt = obj.dt;
            [p_ref, v_ref, a_ref] = get_ref_traj(obj.t);
            x_star = [p_ref; v_ref; 0; 0];
            u_star = 0;
           
            % 1. Linearize DT Dyn about current reference
            A_star = [1, dt, 0, 0;
                 (5*dt*r_g^2*x_star(4)^2*cos(x_star(3))^2)/(7*L^2),  1, (5*dt*r_g*(2*L*g*cos(x_star(3)) + L*r_g*x_star(4)^2*sin(2*x_star(3)) - 2*r_g*x_star(1)*x_star(4)^2*sin(2*x_star(3))))/(14*L^2), -(5*dt*r_g^2*x_star(4)*cos(x_star(3))^2*(L - 2*x_star(1)))/(7*L^2);
                 0,  0,  1, dt;
                 0,  0,  0, 1 - dt/tau]; 
            B_star = [0;0;0;dt*(K/tau)];

            % 2. LQR weights
            Q = diag([1000, 500, 1, 1]);
            R = diag(0.1);

            % 3. Apply LQR
            % [K_lqr,~,~] = dlqr(A_star, B_star, Q, R);
            K_lqr = dlqr_iterative(obj, A_star, B_star, Q, R);
            
            x_est = obj.xm;
            V_servo = u_star - K_lqr*(x_est-x_star);

            % 4. Update parameters
            obj.V_servo = V_servo;
        end

        %% Auxilary functions
        function K = dlqr_iterative(obj, Ad, Bd, Q, R)
    
            tol = 1e-5;
            max_iter = 10000;
            
            % Initialize
            P = Q;
            iter = 0;
            error = 0;
            
            while iter < max_iter
                iter = iter + 1;
                
                P_next = Q + Ad'*P*Ad - Ad'*P*Bd / (R + Bd'*P*Bd) * Bd'*P*Ad;
                
                error = norm(P_next - P, 'fro');
                if error < tol
                    P = P_next;
                    break;
                end
                
                P = P_next;
            end
            
            if iter == max_iter
                fprintf("At %f: DARE iteration did not converge within max iteration. Error: %f\n", obj.t, error);
            end
            
            % Compute the optimal state feedback gain
            K = (R + Bd'*P*Bd) \ (Bd'*P*Ad);
        
        end

    end
    
    methods(Access = public)

        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)        
            % Run setup at start
            % if t==0
            %     setupImpl(obj)
            % end
            % Step controller, observer, etc.
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
    end
    
end