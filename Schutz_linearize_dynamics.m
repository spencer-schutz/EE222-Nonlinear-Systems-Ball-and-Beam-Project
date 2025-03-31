syms x1 x2 x3 x4 u r_g L g K tau real
x = [x1;x2;x3;x4];
dxdt = [x2;
        (5*g/7)*(r_g/L)*sin(x3)-(5/7)*(L/2-x1)*(r_g/L)^2*x4^2*(cos(x3))^2;
        x4;
        -x4/tau + (K/tau)*u];
A = simplify(expand(jacobian(dxdt, x)));
B = simplify(expand(jacobian(dxdt, u)));

