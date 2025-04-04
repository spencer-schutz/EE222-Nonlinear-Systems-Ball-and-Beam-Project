classdef Part1_studentControllerInterface < matlab.System
    properties (Access = private)
        % Time, reference, measurement, input
        t = 0;
        t_prev = -1;
        dt = 0;
        ref = 0;
        meas = 0;
        u = 0;

        % NL CT Dynamics xdot(t)=f(x(t),u(t))
        f_c = @(x, u, r_g, L, g, K, tau) ...
                 [x(2);
                  (5*g/7)*(r_g/L)*sin(x(3))-(5/7)*(L/2-x(1))*(r_g/L)^2*x(4)^2*(cos(x(3)))^2;
                  x(4);
                  -x(4)/tau + (K/tau)*u];
        
        % L CT Dynamics xdot(t)=A(t)x(t)+B(t)u(t)
        A_c = @(x, r_g, L, g, K, tau) ...
               [0, 1, 0, 0;
                (5*r_g^2*x(4)^2*cos(x(3))^2)/(7*L^2), 0, (5*r_g*(2*L*g*cos(x(3)) + L*r_g*x(4)^2*sin(2*x(3)) - 2*r_g*x(1)*x(4)^2*sin(2*x(3))))/(14*L^2), -(5*r_g^2*x(4)*cos(x(3))^2*(L - 2*x(1)))/(7*L^2);
                0, 0, 0, 1;
                0, 0, 0, -1/tau];
        
        % L DT Dynamics x[k]=A[k-1]x[k-1]+B[k-1]u[k-1]
        A_d = @(x, r_g, L, g, K, tau, dt) ...
                [1, dt, 0, 0;
                 (5*dt*r_g^2*x(4)^2*cos(x(3))^2)/(7*L^2),  1, (5*dt*r_g*(2*L*g*cos(x(3)) + L*r_g*x(4)^2*sin(2*x(3)) - 2*r_g*x(1)*x(4)^2*sin(2*x(3))))/(14*L^2), -(5*dt*r_g^2*x(4)*cos(x(3))^2*(L - 2*x(1)))/(7*L^2);
                 0,  0,  1, dt;
                 0,  0,  0, 1 - dt/tau]; 
        B_d = @(K, tau, dt)...
               [0;0;0;dt*(K/tau)];
        L_d = @(dt)...
                dt*eye(4);      

        % Measurement Equation (for EKF)
        H_d = [1 0 0 0; 0 0 1 0];
        M_d = eye(2);

        % EKF variables
        xm = [-.19; 0; 0; 0]; % init with x0
        Pm = diag([0.05, 0, 0.01, 0]); % init with P0

        % PID variables
        theta_d = 0;        
        theta_saturation = 56 * pi / 180;

        % Feedback Linearization variables
        sig_x = 0;
        sig_theta = 0;
        d2traj_ref = [];

        % Error term
        sum_e = 0;
    end

    methods(Access = protected)
        %% Step + Step
        function setupImpl(obj)
            disp("You can use this function for initializaition.");
        end

        function V_servo = stepImpl(obj, t, p_ball, theta)
            % Publish current step info
            [p_ref, v_ref, a_ref] = get_ref_traj(t);
            ref = [p_ref, v_ref, a_ref];
            obj.ref = ref;

            meas = [p_ball; theta];
            obj.meas = meas;
            
            obj.t = t;
            obj.dt = t - obj.t_prev;

            % Choose observer
            %EKF_Carlo(obj)
            %ExtendedHighGain(obj)
            EKF_Spencer(obj)

            % % Publish error
            if obj.t > 3
                e_now = obj.xm(1:2,:) - obj.ref(1:2)';
                obj.sum_e = obj.sum_e + e_now * obj.dt;   
            end

            % Choose controller 
            %FeedbackLinearization(obj)
            LQR(obj)
            %PID(obj)

            % Save output
            V_servo = obj.u;

            % Go to next step
            obj.t_prev = t;
        end

        %% Observers
        function EKF_Spencer(obj)
            % 0. Setup
            r_g = 0.0254;
            L = 0.4255;
            g = 9.81;
            K = 1.5;
            tau = 0.025;

            dt = obj.dt;
            meas = obj.meas;

            Svv = diag([0.05, 0.001, .0175, .00175]);
            Sww = diag([0.05, 0.0175]);

            % 1. Predict state - rollout discretized NL dynamics (forward Euler)
            xm_prev = obj.xm; % from previous iteration 
            u_prev = obj.u; % known
            xp_now = xm_prev + dt*obj.f_c(xm_prev, u_prev, r_g, L, g, K, tau);

            % 2. Predict variance - linearize about last estimate
            Pm_prev = obj.Pm; % from previous iteration 
            A_prev = obj.A_d(xm_prev, r_g, L, g, K, tau, dt);
            L_prev = obj.L_d(dt);
            Pp_now = A_prev*Pm_prev*A_prev' + L_prev*Svv*L_prev';

            % 3. Measurement update
            H_now = obj.H_d;
            M_now = obj.M_d;
            z_now = meas;
            zp_now = H_now*xp_now;
            K_now = Pp_now*H_now'*(H_now*Pp_now*H_now' + M_now*Sww*M_now')^-1;
            xm_now = xp_now + K_now*(z_now - zp_now);
            Pm_now = (eye(4) - K_now*H_now) * Pp_now * (eye(4)-K_now*H_now)' + K_now*M_now*Sww*M_now'*K_now';
                
            % 4. Update Parameters
            obj.xm = xm_now;
            obj.Pm = Pm_now;
        end

        function EKF_Carlo(obj)
            % 0. Setup
            r_g = 0.0254;
            L = 0.4255;
            g = 9.81;
            K = 1.5;
            tau = 0.025;

            x_est = obj.xm;
            u = obj.u;
            dt = obj.dt;
            
            Svv = 0.1*eye(4); % process noise covariance
            Sww = 0.1*eye(2); % measurement noise covariance
             
            % prediction step
            xp = obj.f_c(x_est, u, r_g, L, g, K, tau);
            
            % linearization
            A = obj.A_c(x_est, r_g, L, g, K, tau);
            X = eye(4,4);
            Pp = A*obj.Pm*A' + X*Svv*X';

            % measurement update
            H = obj.H_d;
            Kk = Pp*H'*inv(H*Pp*H' + Sww);
            z = obj.meas;
            obj.xm = xp + Kk*(z - H*xp);
            obj.Pm = (eye(4) - Kk*H)*Pp;
           
        end
        
        function ExtendedHighGain(obj)
            % 0. Setup
            r_g = 0.0254;
            L = 0.4255;
            g = 9.81;
            K = 1.5;
            tau = 0.025;

            a = 5/7*g*r_g/L;
            sat_sig_x = 0.8;
            sat_sig_theta = 0.8;
            
            % Gains
            epsilon = 0.04;
            
            % State parameters
            t_prev = obj.t_prev; 
            dt = obj.dt;
            x_prev = obj.xm;
            sig_x_prev = obj.sig_x; 
            sig_theta_prev = obj.sig_theta;
            u_prev = obj.u;
            meas = obj.meas;
            p_ball = meas(1,1);
            theta = meas(2,1);
                        
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
            obj.xm = x_est;
            obj.sig_x = sig_x; 
            obj.sig_theta = sig_theta;

        end

        %% Controllers
        function FeedbackLinearization(obj)
            % 0. Setup
            r_g = 0.0254;
            L = 0.4255;
            g = 9.81;
            K = 1.5;
            tau = 0.025;

            a = 5/7*g*r_g/L;

            % Fit polynomial to get derivative
            ref = obj.ref;
            p_ref = ref(1);
            dp_ref = ref(2);
            d2p_ref = ref(3);
            d2traj_ref = [obj.d2traj_ref, d2p_ref];
            windowSize = 7;
            polyOrder = 2;
            
            t = obj.t;
            dt = obj.dt;
            
            % Get higher order derivatives of reference
            if t >= 0.1
                p = polyfit(t-dt*windowSize:dt:t, d2traj_ref(end-windowSize:end), polyOrder);
                d3p_ref = polyval(polyder(p), t);
                d4p_ref = polyval(polyder(polyder(p)), t);
            else
                d3p_ref = 0;
                d4p_ref = 0;
            end
            obj.d2traj_ref = d2traj_ref;

            % Approximate Feedback Linearization and DLQR setup
            x_est = obj.xm;
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
            u = 1/g_x * (-f_x + d4p_ref - K_lqr(1)*e1 - K_lqr(2)*e2 - K_lqr(3)*e3 - K_lqr(4)*e4);

            % Update properties
            obj.u = u;
        end

        function LQR(obj)
            % 0. Constants
            r_g = 0.0254;
            L = 0.4255;
            g = 9.81;
            K = 1.5;
            tau = 0.025;

            dt = obj.dt;
            ref = obj.ref;
            x_star = [ref(1); ref(2); 0; 0];
            u_star = 0;
           
            % 1. Linearize DT Dyn about current reference
            A_star = obj.A_d(x_star, r_g, L, g, K, tau, dt);
            B_star = obj.B_d(K, tau, dt);

            % 2. LQR weights
            Q = diag([1000, 500, 1, 1]);
            R = diag(0.1);

            % 3. Apply LQR
            [K_lqr,~,~] = dlqr(A_star, B_star, Q, R);
            x_est = obj.xm;
            u = u_star - K_lqr*(x_est-x_star);
            
            % Augment I term
            if obj.t > 3
                Ke = [150 100];
                u = u - Ke*obj.sum_e;
            end

            % 4. Update parameters
            obj.u = u;
        end

        function PID(obj)
            % 0. Setup
            ref = obj.ref;
            p_ball_ref = ref(1);

            x_est = obj.xm;
            p_est = x_est(1,1);
            v_est = x_est(2,1);

            theta = obj.meas(2,1);
         
            % Decide desired servo angle based on simple proportional feedback.
            k_p = 3;
            k_d = 10;
            obj.theta_d = - k_p * (p_est - p_ball_ref) - k_d*v_est;

            % Input limits
            obj.theta_d = min(obj.theta_d, obj.theta_saturation);
            obj.theta_d = max(obj.theta_d, -obj.theta_saturation);

            % Simple position control to control servo angle to the desired
            % position.
            k_servo = 5;
            u = k_servo * (obj.theta_d - theta);
            
            % Update parameters
            obj.u = u;
        end

    end
    
    methods(Access = public)

        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)        
            % Run setup at start
            if t==0
                setupImpl(obj)
            end
            % Step controller, observer, etc.
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
    end
    
end