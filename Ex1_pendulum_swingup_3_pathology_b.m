% clc;
clear;
rng(10);
%% Parameters
n = 2; % Original state order;
l = 0.15;
m = 0.017;
g = 9.84;
Ts = 0.01;
gamma = 0.001;
% Inverted Pendulum
f_pend = @(t,x,u) [x(2); g/l*sin(x(1)) - gamma*x(2)/(m*l^2) + u/(m*l^2)];

Ntraj = 2000; % Number of trajectories of each identification dataset
Ntraj_val = 1000;
tspan = 0:Ts:0.5; % Time span of each identification dataset

% Final simulation
x0 = [pi;0];
time_vec = 0:Ts:4;

%% Collect Idf data
f_u = @(t, x, traj) 0.05*randn(); % Identification input

% #1
X0 = pi - [0.2*pi*(randn(1, Ntraj)); randn(1, Ntraj)];
[X,Y,U] = collect_trajs(f_pend, f_u, X0, Ntraj, tspan, 0, "ode4");
[Xsep_idf, ~, Usep_idf] = separate_trajs(X,Y,U, Ntraj, tspan(end), Ts);

% #2
X0 = [0.1*pi*(randn(1, Ntraj)); randn(1, Ntraj)];
[X2,Y2,U2] = collect_trajs(f_pend, f_u, X0, Ntraj, tspan, 0, "ode4");
[Xsep2_idf, ~, Usep2_idf] = separate_trajs(X2,Y2,U2, Ntraj, tspan(end), Ts);

% #3
X0 = [pi*(randn(1, Ntraj)); randn(1, Ntraj)];
[X3,Y3,U3] = collect_trajs(f_pend, f_u, X0, Ntraj, tspan, 0, "ode4");
[Xsep3_idf, ~, Usep3_idf] = separate_trajs(X3,Y3,U3, Ntraj, tspan(end), Ts);

%% Gather trajectories for evaluation
X0_eval = [pi-0.1*pi*(randn(1, Ntraj)); randn(1, Ntraj)];
f_u = @(t, x, traj) 0.05*randn();
[X_eval,Y_eval,U_eval] = collect_trajs(f_pend, f_u, X0_eval, Ntraj_val, tspan, 0, "ode4");
[Xsep, Ysep, Usep] = separate_trajs(X_eval,Y_eval,U_eval, Ntraj_val, tspan(end), Ts);

%% Identification
idx = 1;

% 1) Linearized system
A_lin = [0, 1; g/l, -gamma];
B_lin = [0; 1/(m*l^2)];
C_lin = eye(2);
systems{idx}.sys = c2d(ss(A_lin,B_lin,C_lin,0, 0), Ts);
systems{idx}.lift_func = @(X) X;
idx = idx+1;

% LLP No. 1
lift_f1 = @(X) [X; sin(X(1,:))];
[A,B,C]= SystemID_via_EDMD(X,Y,U, lift_f1);
systems{idx}.sys = ss(A,B,C,0, Ts);
systems{idx}.lift_func = lift_f1;
idx = idx+1;

% LLP No. 2
lift_f1 = @(X) [X; sin(X(1,:))];
[A,B,C]= SystemID_via_EDMD(X2,Y2,U2, lift_f1);
systems{idx}.sys = ss(A,B,C,0, Ts);
systems{idx}.lift_func = lift_f1;
idx = idx+1;

% LLP No. 3
lift_f1 = @(X) [X; sin(X(1,:))];
[A,B,C]= SystemID_via_EDMD(X3,Y3,U3, lift_f1);
systems{idx}.sys = ss(A,B,C,0, Ts);
systems{idx}.lift_func = lift_f1;
idx = idx+1;


%% Error computation
Pred_error = cell(1, numel(systems));
Step_errors = cell(1, numel(systems));

for idx_sys = 1:numel(systems)
    % ------ One Step Errors ------
    cost_funcs_step = @(x, u, x_n, A, B, C, lift_f) [
        mean( norm_vec(x_n - C*(A*lift_f(x) + B*u)).^2);
        mean( norm_vec(lift_f(x_n) - (A*lift_f(x) + B*u)).^2);
        ]/Ntraj_val;
    Step_errors{idx_sys} = eval_pred_step_error(X_eval, Y_eval, U_eval, cost_funcs_step, systems{idx_sys});

    % ------ Prediction Error ------
    cost_funcs_pred = @(x, x_eval) mean(norm_vec((x - x_eval)').^2)/Ntraj_val;
    Pred_error{idx_sys} = eval_pred_predict_error(systems{idx_sys}, Xsep, Usep, cost_funcs_pred, tspan);
end

%% Closed-loop for every LLP
XX = cell(1, numel(systems));
UU = cell(1, numel(systems));

for predictor_idx = 1:numel(systems)
    curr_predictor = systems{predictor_idx};

    % Design LQR:
    Q = blkdiag(eye(n), zeros(size(curr_predictor.sys.A,1)-n));
    R = 1000;
    [K,~,~] = dlqr(curr_predictor.sys.A,curr_predictor.sys.B,Q,R,0);
    lift_func = curr_predictor.lift_func;

    % Initialize simulation
    Xlqr = zeros(numel(x0), numel(time_vec));
    Xlqr(:,1) = x0;
    Xprop = Xlqr;
    Usim = zeros(1, numel(time_vec));
    
    % Iterate thought time
    for ii = 2:numel(time_vec)
        Usim(ii-1) = -K*lift_func(Xlqr(:,ii-1));
        Xsim = ode4(@(t,x) f_pend(t,x,Usim(ii-1)), [0, Ts], Xlqr(:,ii-1));
        Xlqr(:,ii) = Xsim(end,:)';
    end
    XX{predictor_idx} = Xlqr;
    UU{predictor_idx} = Usim;
end

%% Generating figures
quant_plot_phase = 10;

gca_fontsize = 12;
MarkerSize_ = 6;
blue = [0 0.447058823529412 0.741176470588235];
dark_red = [0.635294117647059 0.0784313725490196 0.184313725490196];
gold = [0.929411764705882 0.694117647058824 0.125490196078431];
grey = [0.501960784313725 0.501960784313725 0.501960784313725];
green = [0.466666666666667 0.674509803921569 0.188235294117647];
purple = [0.494 0.184 0.556];

set(0,'defaulttextinterpreter','latex');
FigHandle = figure('Position', [100, 100, 600, 450]); % nice proportions for 2-column article

hold on;grid on;box on;
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel

% ====================
subplot(2,3,1);
hold on;grid on;box on;
ylabel("$\omega$ (rad/s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
for ii = 1:quant_plot_phase:Ntraj
    plot(Xsep_idf{ii}(1,1:1:end), Xsep_idf{ii}(2,1:1:end), 'LineStyle',':', 'Color', grey);
end
plot(pi,0, 'MarkerSize',MarkerSize_,'Marker','o','LineWidth',3,'LineStyle','none','Color',gold);
plot(0,0, 'MarkerSize',MarkerSize_,'Marker','o','LineWidth',3,'LineStyle','none','Color',gold);
xlim([-2*pi 2*pi]);
ylim([-20 20]);
title("No. 1", 'fontsize',gca_fontsize, 'Interpreter','LaTex');
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
set(gca,'TickLabelInterpreter','latex');

% ====================
subplot(2,3,2);
hold on;grid on;box on;
xlabel("$\varphi$ (rad)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
for ii = 1:quant_plot_phase:Ntraj
    plot(Xsep2_idf{ii}(1,1:1:end), Xsep2_idf{ii}(2,1:1:end), 'LineStyle',':', 'Color', grey);
end
plot(pi,0, 'MarkerSize',MarkerSize_,'Marker','o','LineWidth',3,'LineStyle','none','Color',gold);
plot(0,0, 'MarkerSize',MarkerSize_,'Marker','o','LineWidth',3,'LineStyle','none','Color',gold);
xlim([-2*pi 2*pi]);
ylim([-20 20]);
title("No. 2", 'fontsize',gca_fontsize, 'Interpreter','LaTex');
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
set(gca,'TickLabelInterpreter','latex');
set(gca,'ytick',[])

% ====================
subplot(2,3,3);
hold on;grid on;box on;
for ii = 1:quant_plot_phase:Ntraj
    plot(Xsep3_idf{ii}(1,1:1:end), Xsep3_idf{ii}(2,1:1:end), 'LineStyle',':', 'Color', grey);
end
plot(pi,0, 'MarkerSize',MarkerSize_,'Marker','o','LineWidth',3,'LineStyle','none','Color',gold);
plot(0,0, 'MarkerSize',MarkerSize_,'Marker','o','LineWidth',3,'LineStyle','none','Color',gold);
xlim([-2*pi 2*pi]);
ylim([-20 20]);
title("No. 3", 'fontsize',gca_fontsize, 'Interpreter','LaTex');
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
set(gca,'TickLabelInterpreter','latex');
set(gca,'ytick',[])

% ====================
subplot(2,3,4:5);
grid on;hold on;box on;
plot(time_vec, XX{1}(1,:), 'LineWidth',2);
plot(time_vec, XX{2}(1,:), 'LineWidth',2, 'LineStyle','--', 'Color', dark_red);
plot(time_vec, XX{3}(1,:), 'LineWidth',2, 'LineStyle',':', 'Color', green);
plot(time_vec, XX{4}(1,:), 'LineWidth',2, 'LineStyle','-.', 'Color', purple);
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
set(gca,'TickLabelInterpreter','latex');
xlabel("Time (s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
ylabel("$\varphi$ (rad)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
xlim([time_vec(1) time_vec(end)]);

h(1) = plot(NaN,NaN,'LineStyle',':', 'Color', grey, 'LineWidth',3);
h(2) = plot(NaN,NaN,'Marker','o','LineWidth',3,'LineStyle','none','Color',gold);
h(3) = plot(NaN,NaN,'LineWidth',3, 'Color',blue);
h(4) = plot(NaN,NaN,'LineWidth',2, 'LineStyle','--', 'Color', dark_red);
h(5) = plot(NaN,NaN,'LineWidth',2, 'LineStyle',':', 'Color', green);
h(6) = plot(NaN,NaN,'LineWidth',2, 'LineStyle','-.', 'Color', purple);
legend(h, "Training Trajectories", "System's equilibria", "Linear Approx.", ...
    "Predictor No. 1", "Predictor No. 2", "Predictor No. 3",'Interpreter','latex');





