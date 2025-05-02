syms x1 x2 x3 x4 u r_g L g K tau v1 v2 v3 v4 w1 w2 dt real
x_t = [x1;x2;x3;x4];
w = [w1;w2];
v = [v1;v2;v3;v4];

% CT NL Dyn
dxdt = [x2 + v1;
        (5*g/7)*(r_g/L)*sin(x3) - (5/7) * (L/2-x1) * (r_g/L)^2 * x4^2 * (cos(x3))^2 + v2;
        x4 + v3;
        -x4/tau + (K/tau)*u + v4];
% DT NL Dyn
x_k = x_t + dt*dxdt;

% DT Linearized Dyn
A = simplify(expand(jacobian(x_k, x_t)));
B = simplify(expand(jacobian(x_k, u)));
L = simplify(expand(jacobian(x_k, v)));

% CT/DT NL Meas
z = [1 0 0 0; 0 0 1 0]*x_t + [w1;w2];

% CT/DT Linearized Meas
H = simplify(expand(jacobian(z, x_t)));
M = simplify(expand(jacobian(z, w)));

% CTRB/OBSV of DT Lin Dyn
P = [B A*B A*A*B A*A*A*B];
Q = [H; H*A; H*A*A; H*A*A*A];
rank(P)
rank(Q)

% % Graveyard
% A_est = [1, dt, 0, 0;
%         (5*dt*r_g^2*xm_now(4)^2*cos(xm_now(3))^2)/(7*L^2),  1, (5*dt*r_g*(2*L*g*cos(xm_now(3)) + L*r_g*xm_now(4)^2*sin(2*xm_now(3)) - 2*r_g*xm_now(1)*xm_now(4)^2*sin(2*xm_now(3))))/(14*L^2), -(5*dt*r_g^2*xm_now(4)*cos(xm_now(3))^2*(L - 2*xm_now(1)))/(7*L^2);
%          0,  0,  1, dt;
%          0,  0,  0, 1 - dt/tau]; % linearize DT NL dyn about xm_prev
% B_est = [0;0;0;dt*(K/tau)];