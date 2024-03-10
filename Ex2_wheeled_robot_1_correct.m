% Author: Loi Do, 2024
% doloi@fel.cvut.cz

clear;
clc;
rng(10);

addpath("theHood/");
addpath("theHood/ode_solvers/");
addpath(genpath("theHood/MPC_functions/"));
%% Parameters
Ts = 0.01; % Sampling period of the resulting linerized system
ulim = 1; % Maximal amplitude of the input (torque)

% Parameters of the input for identification/validation trajectories
Ntraj = 1000; % Number of identification trajectories


%% Collect Idf data

% Linearize system
x0 = zeros(6,1); % Operation point
u0 = zeros(2,1);
[A_lin,B_lin] = linearizeSys(@(t,x,u)f_3DwheeledRob(t,x,u), x0, u0, Ts);
[A_lin_reduc,B_lin_reduc] = remove_states(A_lin, B_lin); % Remove states: distance, yaw
C_lin = eye(size(A_lin_reduc,1));

% LQR weight matrices
Q = diag([20,0.237,0.119,20]);
R = 10*eye(2);
[K,S,e] = dlqr(A_lin_reduc,B_lin_reduc,Q,R);

f_u = @(t,x, ref, traj)  bound_val(-K*(x - ref) + 0.1*rndBx(2,1), -ulim, ulim); 
tspan = 0:Ts:0.5;
t_size = numel(tspan);
X0 = 0.25*rndBx(4,Ntraj);

% References
dot_S_Ref = (2*(rand(Ntraj,1)-0.5)*tspan)';
dot_Phi_Ref = zeros(1, Ntraj*t_size);
dot_Psi_Ref = ( 2*(rand(Ntraj,1)-0.5)*tspan)';
Ref = [dot_S_Ref(:)'; dot_Phi_Ref; dot_Psi_Ref(:)'; zeros(1, Ntraj*t_size)];

f_sys = @(t,x,u) f_3DwheeledRob_reduc(t,x,u);
[X,Y,U] = collect_trajs_ref(f_sys, f_u, X0, Ref, Ntraj, tspan, 0, "ode4");
%% Gather trajectories for evaluation
Ntraj_eval = 1000;
[X_eval,Y_eval,U_eval] = collect_trajs_ref(f_sys, f_u, X0, Ref, Ntraj_eval, tspan, 0, "ode4");
[Xsep, Ysep, Usep] = separate_trajs(X_eval,Y_eval,U_eval, Ntraj_eval, tspan(end), Ts);

%% Identification
idx = 1;

% 1) Linearized system
systems{idx}.sys = ss(A_lin_reduc,B_lin_reduc,C_lin,0);
systems{idx}.lift_func = @(X) X;
idx = idx+1;

% 2) Lifting with sinus
lift_f1 = @(X) [X; sin(X(4,:));];
[A,B,C]= SystemID_via_EDMD(X,Y,U, lift_f1);
systems{idx}.sys = ss(A,B,C,0, Ts);
systems{idx}.lift_func = lift_f1;
idx = idx+1;

%% Error computation
Pred_error = cell(1, numel(systems));
Step_errors = cell(1, numel(systems));

% Eval
for idx_sys = 1:numel(systems)
    % ------ Prediction Errors ------
    cost_funcs_pred = @(x, x_eval) [
        mean(norm_vec((x - x_eval)').^2);
        ]/Ntraj_eval;
    Pred_error{idx_sys} = eval_pred_predict_error(systems{idx_sys}, Xsep, Usep, cost_funcs_pred, tspan);

    % ------ Step Errors ------
    cost_funcs_step = @(x, u, x_n, A, B, C, lift_f) [
        mean( norm_vec(x_n - C*(A*lift_f(x) + B*u)).^2);
        mean( norm_vec(lift_f(x_n) - (A*lift_f(x) + B*u)).^2);
        ]/Ntraj_eval;
    Step_errors{idx_sys} = eval_pred_step_error(X_eval, Y_eval, U_eval, cost_funcs_step, systems{idx_sys});
end
Step_errors_mat = cell2mat(Step_errors);
Pred_error_mat = cell2mat(Pred_error);

%% Closed-loop simulation
XX = cell(1, numel(XX0));
UU = cell(1, numel(XX0));

% Reference for simulation
time_vec = 0:Ts:10;
t_size = numel(time_vec);
dot_S_Ref_Comp = -0.8*((square(1*time_vec))-1);
dot_Psi_Ref_Comp = -3*(square(3*time_vec));
dot_Phi_Ref_Comp = zeros(1, t_size); % Keeping at zero ideally
Ref_Comp = [dot_S_Ref_Comp; dot_Phi_Ref_Comp; dot_Psi_Ref_Comp; zeros(1, t_size)];

x0 = zeros(4,1);
for predictor_idx = 1:numel(systems)
    best_pred = systems{predictor_idx};

    % Design LQR:
    Q_lift = blkdiag(Q, zeros(size(best_pred.sys.A,1)-size(Q,1)));

    Alift = best_pred.sys.A;
    Blift = best_pred.sys.B;
    Clift = best_pred.sys.C;
    lift_func = best_pred.lift_func;
    n =  size(Alift,1);

    % KMPC
    Npred = 50;
    xlift_min = -1e10*ones(n,1); % Lower bounds on lifted states
    xlift_max = 1e10*ones(n,1); % Upper bounds on lifted states
    umin_mpc = -ulim*ones(2,1);
    umax_mpc = ulim*ones(2,1);
    [~,~,kmpc]  = osqp_MPC_controller(Alift,Blift,Clift,Q,R,Q,Npred,umin_mpc, umax_mpc, xlift_min, xlift_max);

    Xsim = zeros(numel(x0), t_size);
    Xsim(:,1) = x0;
    Usim = zeros(2, t_size);

    for ii = 2:t_size
        z = lift_func(Xsim(:,ii-1));
        try
            [u, ~] = kmpc(z,Ref_Comp(:,ii-1)); % get the control input
        catch
            break;
        end
        Usim(:,ii-1) = u(1:2);
        Xsim_curr = ode4(@(t,x) f_3DwheeledRob_reduc(t,x,Usim(:,ii-1)), [0, Ts], Xsim(:,ii-1));
        Xsim(:,ii) = Xsim_curr(end,:)';
    end
    XX{predictor_idx} = Xsim;
    UU{predictor_idx} = Usim;
end

%% Plot figures

gca_fontsize = 12;
set(0,'defaulttextinterpreter','latex');
FigHandle = figure('Position', [100, 100, 600, 450]); % nice proportions for 2-column article

% Cut-off for trajectories
idx_end = [720, numel(time_vec)];
for ii = 1:numel(systems)

    subplot(4,1,1);
    hold on;
    grid on;
    box on;
    plot(time_vec(1:idx_end(ii)),XX{ii}(4,1:idx_end(ii)),'Color',color_p(ii), 'LineWidth',2);
    ylabel("$\varphi$ (rad)",'FontSize',gca_fontsize, 'Interpreter','LaTex');

    subplot(4,1,2);
    hold on;
    grid on;
    box on;
    plot(time_vec(1:idx_end(ii)),XX{ii}(1,1:idx_end(ii)), 'Color',color_p(ii), 'LineWidth',2);
    ylabel("$\dot{s}$ (m/s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');

    subplot(4,1,3);
    hold on;
    grid on;
    box on;
    plot(time_vec(1:idx_end(ii)),XX{ii}(2,1:idx_end(ii)),'Color',color_p(ii), 'LineWidth',2); 
    ylabel("$\dot{\varphi}$ (rad/s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');

    subplot(4,1,4);
    hold on;
    grid on;
    box on;
    plot(time_vec(1:idx_end(ii)),XX{ii}(3,1:idx_end(ii)),'Color',color_p(ii), 'LineWidth',2);
    ylabel("$\dot{\chi}$ (rad/s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
    xlabel("Time (s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');

end

% Plot end markers
idx_data = [4,1,2,3];
for ii = 1:4
    subplot(4,1,ii);
    plot(time_vec(idx_end(1)),XX{1}(idx_data(ii),idx_end(1)),'Color',color_p(1), 'LineWidth',5, 'Marker','o', 'MarkerSize',3);
end


subplot(4,1,1);
plot(time_vec,Ref_Comp(4,:), 'LineWidth',2, "LineStyle","--","Color", "black");
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',gca_fontsize);
ylim([-1.5 0.8]);

subplot(4,1,2);
plot(time_vec,Ref_Comp(1,:), 'LineWidth',2, "LineStyle","--","Color", "black");
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',gca_fontsize);
ylim([-1 2]);

subplot(4,1,3);
plot(time_vec,Ref_Comp(2,:), 'LineWidth',2, "LineStyle","--","Color", "black");
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',gca_fontsize);
ylim([-5 2]);

subplot(4,1,4);
plot(time_vec,Ref_Comp(3,:), 'LineWidth',2, "LineStyle","--","Color", "black");
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',gca_fontsize);
ylim([-6 5]);

h(1) = plot(NaN,NaN,'Color',color_p(1), 'LineWidth',2);
h(2) = plot(NaN,NaN,'Color',color_p(2), 'LineWidth',2);
h(3) = plot(NaN,NaN,'Color',"black", 'LineWidth',2, "LineStyle","--");

legend(h, "Local at $x_0$", "Lifted Predictor", "Reference", 'Position',[0.737421888487879 0.016881238209532 0.253133667067676 0.126311115434435],...
    'Interpreter','latex');

xlim([0, 10]);







