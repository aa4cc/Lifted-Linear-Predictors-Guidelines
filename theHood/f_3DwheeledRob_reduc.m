function dxdt = f_3DwheeledRob_reduc(t,x, u)
%F_3DWHEELEDROB ODE of 3D wheeled robot.
% u(t,x) should be a anonymous function with two outputs (motors' torques)
% and inputs: 't': current sim time, 'x': current state

% System's parameters
b = 1*10^-2;
J = 7.35*10^-4;
K = 3.9*10^-4;
m_w = 3*10^-1;
r = 8*10^-2;
w = 2.9*10^-1;
l = 2.907*10^-1;
I_px = 1.5625*10^-2;
I_py = 1.18625*10^-2;
I_pz = 1.18625*10^-2;
m_p = 4;
g = 9.81;

lims = [1; 1];

% State variables
dots = x(1);    % Traveled speed
dotphi = x(2);  % Tilt angular speed
dotpsi = x(3);  % Yaw angular speed
% s = 0;   % Distance traveled
phi = x(4); % Tilt angle
% psi = 0; % Yaw angle


m11 = m_p + 2*m_w + 2*J/r^2;
m12 = m_p*l*cos(phi);
m21 = m12;
m22 = I_py + m_p*l^2;
m33 = I_pz + 2*K + (m_w + J/r^2)*w^2/2 - (I_pz - I_px - m_p*l^2)*sin(phi)^2;

c12 = -m_p*l*dotphi*sin(phi);
c13 = m_p*l*dotpsi*sin(phi);
c23 = (I_pz - I_px - m_p*l^2)*dotpsi*sin(phi)*cos(phi);
c31 = m_p*l*dotpsi*sin(phi);
c32 = -c23;
c33 = -(I_pz - I_px - m_p*l^2)*dotphi*sin(phi)*cos(phi);

d11 = 2*b/r^2;
d12 = -2*b/r;
d21 = d12;
d22 = 2*b;
d33 = w^2*b/(2*r^2);


M = [m11, m12, 0; 
     m21, m22, 0;
     0, 0, m33];
C = [0, c12, c13; 
     0, 0, c23;
     c31, c32, c33];
D = [d11, d12, 0;
     d21, d22, 0;
     0, 0, d33];
B = [-1/r, -1/r;
    1, 1;
    w/(2*r), -w/(2*r)];
G = [0; -m_p*l*g*sin(phi); 0];



torque_in = max(min(u, lims), -lims); 


dotq = [dots; dotphi; dotpsi];
dotdotq = M\(B*torque_in - (C + D)*dotq - G);

dxdt = [
    dotdotq;
    dotq(2)
];


end

