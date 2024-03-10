function [A,B] = linearizeSys(f_sys, x0, u0, Ts)
%linearizeSys : Linearize the system. Tailored for 3D (Sk8o) robot but
%might be usable for other systems
%
% f_sys(t,x,u): function to be linearized
% Ts: Sampling period of the resulting system, Ts = 0: for continuous
% x0: linearization point
% u0: linearization input
% system.
% Return: A,B matrices of the resulting system

% Author: Loi Do, 2024

x_sym = sym('x',[numel(x0),1]);
u_sym = sym('u',[numel(u0),1]);

f = f_sys(0,x_sym, u_sym);
[A,B] = jacobians(f,x_sym,u_sym,x0,u0);

A_c = double(A);
B_c = double(B);
C_c = eye(6);
sys_c = ss(A_c, B_c, C_c, 0);

if Ts == 0
    A = A_c;
    B = B_c;
else
    sys_d = c2d(sys_c, Ts);
    A = sys_d.A;
    B = sys_d.B;
end

end