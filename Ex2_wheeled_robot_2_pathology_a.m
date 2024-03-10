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
C_lin = eye(size(A_lin,1));

% LQR weight matrices
Q = diag([20,0.237,0.119,30,20,1]);
R = 10*eye(2);
[K,S,e] = dlqr(A_lin,B_lin,Q,R);

f_u = @(t,x, ref, traj) bound_val(-K*(x - ref) + 0.1*rndBx(2,1), -ulim, ulim);
tspan = 0:Ts:0.5;
X0 = 0.25*rndBx(6,Ntraj);

% References
dot_S_Ref = (2*(rand(Ntraj,1)-0.5)*tspan)';
dot_Phi_Ref = zeros(1, Ntraj*numel(tspan));
dot_Psi_Ref = ( 2*(rand(Ntraj,1)-0.5)*tspan)';
Ref_V = [dot_S_Ref(:)'; dot_Phi_Ref; dot_Psi_Ref(:)'];
Ref = [Ref_V; Ts*cumsum(Ref_V)];

f_sys = @(t,x,u) f_3DwheeledRob(t,x,u);
[X,Y,U] = collect_trajs_ref(f_sys, f_u, X0, Ref, Ntraj, tspan, 0, "ode4");
%% Gather trajectories for evaluation
Ntraj_eval = 1000;
[X_eval,Y_eval,U_eval] = collect_trajs_ref(f_sys, f_u, X0, Ref, Ntraj_eval, tspan, 0, "ode4");
[Xsep, Ysep, Usep] = separate_trajs(X_eval,Y_eval,U_eval, Ntraj_eval, tspan(end), Ts);


%% Identification
clear systems
idx = 1;
lift_f1 = @(X) [X; sin(X(4,:))];    % LLP with Sine
[A,B,C]= SystemID_via_EDMD(X,Y,U, lift_f1);
systems{idx}.sys = ss(A,B,C,0, Ts);
systems{idx}.lift_func = lift_f1;

%% Closed-loop simulation


% References
time_vec = 0:Ts:10;
dot_S_Ref_Comp = -0.8*((square(1*time_vec))-1);
dot_Psi_Ref_Comp = -3*(square(3*time_vec));
dot_Phi_Ref_Comp = zeros(1, numel(time_vec)); % Keeping at zero ideally

s0 = 0;
yaw0 = 0;

Ref_Comp{1} = [dot_S_Ref_Comp; dot_Phi_Ref_Comp; dot_Psi_Ref_Comp; Ts*cumsum(dot_S_Ref_Comp)+s0;zeros(1, numel(time_vec));yaw0+Ts*cumsum(dot_Psi_Ref_Comp)];
XX0{1} = zeros(6,1);
XX0{1}(4) = s0;
XX0{1}(6) = yaw0;

s0 = 100;
yaw0 = 0;
XX0{2} = zeros(6,1);
XX0{2}(4) = s0;
XX0{2}(6) = yaw0;
Ref_Comp{2} = [dot_S_Ref_Comp; dot_Phi_Ref_Comp; dot_Psi_Ref_Comp; Ts*cumsum(dot_S_Ref_Comp)+s0;zeros(1, numel(time_vec));yaw0+Ts*cumsum(dot_Psi_Ref_Comp)];


XX = cell(1, numel(XX0));
UU = cell(1, numel(XX0));

predictor = systems{1};

Alift = predictor.sys.A;
Blift = predictor.sys.B;
Clift = predictor.sys.C;
lift_func = predictor.lift_func;
n =  size(Alift,1);

% KMPC
Npred = 50;
xlift_min = -1e10*ones(n,1); % Lower bounds on lifted states
xlift_max = 1e10*ones(n,1); % Upper bounds on lifted states
umin_mpc = -ulim*ones(2,1);
umax_mpc = ulim*ones(2,1);
[~,~,kmpc]  = osqp_MPC_controller(Alift,Blift,Clift,Q,R,Q,Npred,umin_mpc, umax_mpc, xlift_min, xlift_max);

for idx = 1:numel(XX0)
    Xsim = zeros(numel(XX0{idx}), numel(time_vec));
    Xsim(:,1) = XX0{idx};
    Usim = zeros(2, numel(time_vec));

    for ii = 2:numel(time_vec)
        z = lift_func(Xsim(:,ii-1));
        try
            [u, ~] = kmpc(z,Ref_Comp{idx}(:,ii-1)); % get the control input
        catch
            break;
        end
        Usim(:,ii-1) = u(1:2);

        Xsim_curr = ode4(@(t,x) f_3DwheeledRob(t,x,Usim(:,ii-1)), [0, Ts], Xsim(:,ii-1));
        Xsim(:,ii) = Xsim_curr(end,:)';
    end
    XX{idx} = Xsim;
    UU{idx} = Usim;
end
%% Plot figures

gca_fontsize = 12;
set(0,'defaulttextinterpreter','latex');
FigHandle = figure('Position', [100, 100, 600, 450]); % nice proportions for 2-column article

idx_end = [numel(time_vec), 127];
for ii = 1:numel(XX0)

    subplot(2,3,1);
    hold on;grid on;box on;
    plot(time_vec(1:idx_end(ii)),XX{ii}(4,1:idx_end(ii)),'Color',color_p(ii), 'LineWidth',2);
    ylabel("$q$",'FontSize',gca_fontsize, 'Interpreter','LaTex');

    subplot(2,3,2);
    hold on;grid on;box on;
    plot(time_vec(1:idx_end(ii)),XX{ii}(5,1:idx_end(ii)), 'Color',color_p(ii), 'LineWidth',2);

    subplot(2,3,3);
    hold on;grid on;box on;
    plot(time_vec(1:idx_end(ii)),XX{ii}(6,1:idx_end(ii)),'Color',color_p(ii), 'LineWidth',2);

    subplot(2,3,4);
    hold on;grid on;box on;
    plot(time_vec(1:idx_end(ii)),XX{ii}(1,1:idx_end(ii)),'Color',color_p(ii), 'LineWidth',2);
    ylabel("$\dot{q}$",'FontSize',gca_fontsize, 'Interpreter','LaTex');

    subplot(2,3,5);
    hold on;grid on;box on;
    plot(time_vec(1:idx_end(ii)),XX{ii}(2,1:idx_end(ii)),'Color',color_p(ii), 'LineWidth',2);
    xlabel("Time (s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');

    subplot(2,3,6);
    hold on;grid on;box on;
    plot(time_vec(1:idx_end(ii)),XX{ii}(3,1:idx_end(ii)),'Color',color_p(ii), 'LineWidth',2);
end

% Plot end markers
idx_data = [4,5,6,1,2,3];
for ii = 1:6
    subplot(2,3,ii);
    plot(time_vec(idx_end(2)),XX{2}(idx_data(ii),idx_end(2)),'Color',color_p(2), 'LineWidth',5, 'Marker','o', 'MarkerSize',3);
end

subplot(2,3,1);
plot(time_vec,Ref_Comp{1}(4,:), 'LineWidth',2, "LineStyle","--","Color", "black");
plot(time_vec,Ref_Comp{2}(4,:), 'LineWidth',2, "LineStyle","--","Color", "black");
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
title("Distance (m)");
xlim([0, time_vec(end)]);
ylim([0, 110]);

subplot(2,3,2);
plot(time_vec,Ref_Comp{1}(5,:), 'LineWidth',2, "LineStyle","--","Color", "black");
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
title("Tilt (rad)");
xlim([0, time_vec(end)]);

subplot(2,3,3);
plot(time_vec,Ref_Comp{1}(6,:), 'LineWidth',2, "LineStyle","--","Color", "black");
plot(time_vec,Ref_Comp{2}(6,:), 'LineWidth',2, "LineStyle","--","Color", "black");
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
title("Yaw (rad)");
xlim([0, time_vec(end)]);
ylim([-6 0.5]);

subplot(2,3,4);
plot(time_vec,Ref_Comp{1}(1,:), 'LineWidth',2, "LineStyle","--","Color", "black");
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
title("Speed (m/s)");
xlim([0, time_vec(end)]);
ylim([-4 2]);

subplot(2,3,5);
plot(time_vec,Ref_Comp{1}(2,:), 'LineWidth',2, "LineStyle","--","Color", "black");
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
title("Tilt Rate (rad/s)");
xlim([0, time_vec(end)]);
ylim([-2 8.5]);

subplot(2,3,6);
plot(time_vec,Ref_Comp{1}(3,:), 'LineWidth',2, "LineStyle","--","Color", "black");
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
title("Yaw Rate (rad/s)");
xlim([0, time_vec(end)]);
ylim([-13 4]);

h(1) = plot(NaN,NaN,'Color',color_p(1), 'LineWidth',2);
h(2) = plot(NaN,NaN,'Color',color_p(2), 'LineWidth',2);
h(3) = plot(NaN,NaN,'Color',"black", 'LineWidth',2, "LineStyle","--");
legend(h, "$s_0 = 0\,$m", "$s_0 = 100\,$m", "Reference", 'Position',[0.737421888487879 0.016881238209532 0.253133667067676 0.126311115434435],...
    'Interpreter','latex');



