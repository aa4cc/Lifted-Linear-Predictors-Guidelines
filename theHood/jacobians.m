function [A,B] = jacobians(xdot,x,u,x0,u0)
%JACOBIANS Linear approximation matrices of a non-linear dynamic model.
%   Convert from xdot = f(x,u) to xdot = A*x + B*u at x = x0, u = u0.
%
%   Syntax:
%      [A,B] = jacobians(f,x,u,x0,u0)
%
%   Example:
%      % Linearization of the simple pendulum model
%      syms theta omega tau m l g                 % Symbols to use
%      u = [tau];                                 % Input vector
%      x = [theta;omega];                         % State vector
%      f = [omega;1/(m*l^2)*tau-g/l*sin(theta)];  % Transition function
%      x0 = [0;0];                                % Reference state
%      u0 = [0];                                  % Reference input
%      [A,B] = jacobians(f,x,u,x0,u0)
%
%   Requires:
%      Symbolic Math Toolbox
%
%   Author:
%      Ildeberto de los Santos Ruiz
%      idelossantos@ittg.edu.mx
%
%   See also EQUILIBRIUM.

A = sym(zeros(numel(x),numel(x)));
B = sym(zeros(numel(x),numel(u)));
for i = 1:numel(x)
    for j = 1:numel(x)
        A(i,j) = diff(xdot(i),x(j));
    end
end
for i = 1:numel(x)
    for j = 1:numel(u)
        B(i,j) = diff(xdot(i),u(j));
    end
end
if nargin > 3
    A = subs(A,[x;u],[x0;u0]);
    B = subs(B,[x;u],[x0;u0]);
end