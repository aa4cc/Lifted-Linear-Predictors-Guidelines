function [Xsep, Ysep, Usep] = separate_trajs(X,Y,U, Ntraj, Tsim, Ts)
% Separate trajectories (X,Y,U) of EDMD algorithm into cells
Xsep = cell(Ntraj, 1);
Ysep = cell(Ntraj, 1);
Usep = cell(Ntraj, 1);

Nstep = Tsim/Ts;

for ii = 1:Ntraj
    idx_start = Nstep*(ii-1) + 1;
    idx_end = Nstep*(ii);
    
    Xsep{ii} = X(:, idx_start:idx_end);
    Ysep{ii} = Y(:, idx_start:idx_end);
    Usep{ii} = U(:, idx_start:idx_end);   
end

end

