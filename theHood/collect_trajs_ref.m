function [X,Y,U] = collect_trajs_ref(f_sys,f_input, X0, Ref, Ntraj, tspan, verbose, solver)
%COLLECT_TRAJS Summary of this function goes here
%   Detailed explanation goes here

Nsim = numel(tspan)-1; 
n = numel(X0(:,1));
p = numel(f_input(0,Ref(:,1), X0(:,1),1));

X = zeros(n, Ntraj*Nsim);
Y = zeros(n, Ntraj*Nsim);
U = zeros(p, Ntraj*Nsim);

X_traj = zeros(n, Nsim);
U_traj = zeros(p, Nsim);

for traj = 1:Ntraj
    if verbose == 1
        disp("Traj: " + num2str(traj))
    end


    idx_start_R = (traj-1)*numel(tspan) + 1;
    idx_end_R = (traj)*numel(tspan);

    Ref_curr = Ref(:,idx_start_R:idx_end_R);
    x_curr = X0(:,traj);
    X_traj(:, 1) = x_curr;
    for t_idx = 2:numel(tspan)
        u_curr = f_input(t_idx,x_curr,Ref_curr(:,t_idx-1),traj);
        tspan_curr = tspan(t_idx-1: t_idx);
        
        if solver == "ode45"
            [~, Xsim] = ode45(@(t,x) f_sys(t,x, u_curr), tspan_curr, x_curr');
        elseif solver == "ode4"
            Xsim = ode4(@(t,x) f_sys(t,x, u_curr), tspan_curr, x_curr);
        elseif solver == "ode5"
            Xsim = ode5(@(t,x) f_sys(t,x, u_curr), tspan_curr, x_curr);
        else
            [~, Xsim] = ode45(@(t,x) f_sys(t,x, u_curr), tspan_curr, x_curr');
        end
        x_curr = Xsim(end, :)';

        % Save 
        X_traj(:,t_idx) = x_curr;
        U_traj(:,t_idx-1) = u_curr;
    end

    idx_start = (traj-1)*Nsim + 1;
    idx_end = (traj)*Nsim;

    U(:,idx_start:idx_end) = U_traj(:,1:end);
    X(:,idx_start:idx_end) = X_traj(:, 1:end-1);
    Y(:,idx_start:idx_end) = X_traj(:, 2:end);
end



end

