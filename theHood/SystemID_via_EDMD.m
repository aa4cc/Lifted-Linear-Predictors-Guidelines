function [A,B,C]= SystemID_via_EDMD(X,Y,U, BuildKoopmanState)
% Extended library with: sin(x), cos(x), v*sin(x), v*cos(x)


%% Extended Dynamic Mode Decomposition


% Stack observables
Xp = BuildKoopmanState(X);
Yp = BuildKoopmanState(Y);
% Xp = X;
% Yp = Y;
Up = U;

% compute the Koopman linear system 
% tic
% disp('Running EDMD ...')
W = [Xp;Up]*[Xp;Up]';
V = Yp*[Xp;Up]';
M = V*pinv(W);

Nlift = size(Xp,1);
A = M(:,1:Nlift);
B = M(:,Nlift+1:end);

nm= size(X,1);
C = zeros(nm,size(A,1));
C(:,1:nm)=eye(nm);
% toc

%% Build Koopman State from the current+ state 'x' for Closed-loop simulations
% BuildKoopmanState = @(x) [x; 
%     sine_X(x); cosine_X(x);
%     v_sine_X(x); v_cosine_X(x)];

end
