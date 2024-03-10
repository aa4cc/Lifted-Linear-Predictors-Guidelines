addpath("theHood/");
addpath("theHood/ode_solvers/")

clear;
rng(10);
%% Parameters
n = 2;          % Order of the original system
l = 0.15;       % Pendulum length
m = 0.017;      % Pendulum mass
g = 9.84;       % Gravity
Ts = 0.01;      % Sampling time
gamma = 0.001;  % Damping coefficient
% Inverted Pendulum
f_pend = @(t,x,u) [x(2); g/l*sin(x(1)) - gamma*x(2)/(m*l^2) + u/(m*l^2)];

% Parameters of the input for identification/validation trajectories
Tsim = 0.5; % Simulation time of identification and evaluation data
Ntraj = 2000;       % Number of identification trajectories
Ntraj_val = 100;   % Number of validation trajectories
AMP = 0.01;
X0_amp = 0.3*pi;
tspan = 0:Ts:Tsim;

% Final simulation
x0 = [pi;0];
time_vec = 0:Ts:3;

%% Collect Identification data
X0 = [pi-X0_amp*(randn(1, Ntraj)); randn(1, Ntraj)];
u_amp = AMP*rand(1, Ntraj);
f_u = @(t, x, traj) 0.001*randn(1) + u_amp(traj)*(-x(1,:));

[X,Y,U] = collect_trajs(f_pend, f_u, X0, Ntraj, tspan, 0, "ode4");
[Xsep_idf, ~, Usep_idf] = separate_trajs(X,Y,U, Ntraj, tspan(end), Ts);

%% Gather trajectories for validation
X0_eval = [pi-X0_amp*(randn(1, Ntraj_val)); randn(1, Ntraj_val)];
u_amp = AMP*rand(1, Ntraj_val);
f_u = @(t, x, traj) 0.001*randn(1) + u_amp(traj)*(-x(1,:));

[X_eval,Y_eval,U_eval] = collect_trajs(f_pend, f_u, X0_eval, Ntraj_val, tspan, 0, "ode4");
[Xsep, Ysep, Usep] = separate_trajs(X_eval,Y_eval,U_eval, Ntraj_val, Tsim, Ts);

%% Identification
% 1) Linearized system
A_lin = [0, 1; g/l, -gamma];
B_lin = [0; 1/(m*l^2)];
C_lin = eye(2);
systems{1}.sys = c2d(ss(A_lin,B_lin,C_lin,0, 0), Ts);
systems{1}.lift_func = @(X) X;

% 2) Just sinus
lift_f1 = @(X) [X; sin(X(1,:))];
[A,B,C]= SystemID_via_EDMD(X,Y,U, lift_f1);
systems{2}.sys = ss(A,B,C,0, Ts);
systems{2}.lift_func = lift_f1;

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

%% Closed-loop for every Predictor
XX = cell(1, numel(systems));
UU = cell(1, numel(systems));

% For each predictor
for predictor_idx = 1:numel(systems)
    curr_predictor = systems{predictor_idx};
    Alift = curr_predictor.sys.A;
    Blift = curr_predictor.sys.B;
    lift_func = curr_predictor.lift_func;
    
    % Design LQR:
    Q = blkdiag(eye(n), zeros(size(Alift,1)-n));
    R = 1000; % Cost matrix for input
    [K,~,~] = dlqr(Alift,Blift,Q,R,0);
    
    % Initialize simulation
    Xlqr = zeros(numel(x0), numel(time_vec));
    Xlqr(:,1) = x0;
    Usim = zeros(1, numel(time_vec));

    for ii = 2:numel(time_vec)
        z = lift_func(Xlqr(:,ii-1));
        Usim(ii-1) = -K*z;
        Xsim = ode4(@(t,x) f_pend(t,x,Usim(ii-1)), [0, Ts], Xlqr(:,ii-1));
        Xlqr(:,ii) = Xsim(end,:)';
    end
    XX{predictor_idx} = Xlqr;
    UU{predictor_idx} = Usim;
end


%% Generating figures
% Parameters of figures
gca_fontsize = 12;
set(0,'defaulttextinterpreter','latex');
% Defining colors
blue = [0 0.447058823529412 0.741176470588235];
dark_red = [0.635294117647059 0.0784313725490196 0.184313725490196];
gold = [0.929411764705882 0.694117647058824 0.125490196078431];
grey = [0.501960784313725 0.501960784313725 0.501960784313725];

FigHandle = figure('Position', [100, 100, 600, 450]);

subplot(3,2,[1 3 5]);
hold on;grid on;box on;
for ii = 1:1:Ntraj
    plot(Xsep_idf{ii}(1,1:1:end), Xsep_idf{ii}(2,1:1:end), 'LineStyle',':', 'Color', grey);
end
plot(XX{1}(1,:), XX{1}(2,:),'LineWidth',3, 'Color',blue);
plot(XX{2}(1,:), XX{2}(2,:),'LineWidth',3,'Color',dark_red, 'LineStyle','--');
xlabel("$\varphi$ (rad)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
ylabel("$\omega$ (rad/s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');

plot(pi,0, 'MarkerSize',8,'Marker','o','LineWidth',3,'LineStyle','none','Color',gold);
plot(0,0, 'MarkerSize',8,'Marker','o','LineWidth',3,'LineStyle','none','Color',gold);
xlim([-pi/2, 7]);
ylim([-25, 15]);
set(gca,'fontsize',gca_fontsize);
set(gca,'TickLabelInterpreter','latex');

subplot(3,2,2);
grid on;hold on;box on;
plot(time_vec, XX{1}(1,:), 'LineWidth',2);
plot(time_vec, XX{2}(1,:), 'LineWidth',2, 'LineStyle','--', 'Color', dark_red);
set(gca,'fontsize',gca_fontsize);
set(gca,'TickLabelInterpreter','latex');
ylabel("$\varphi$ (rad)",'FontSize',gca_fontsize, 'Interpreter','LaTex');

subplot(3,2,4);
grid on;hold on;box on;
plot(time_vec, XX{1}(2,:), 'LineWidth',2);
plot(time_vec, XX{2}(2,:), 'LineWidth',2, 'LineStyle','--', 'Color', dark_red);
set(gca,'fontsize',gca_fontsize);
set(gca,'TickLabelInterpreter','latex');
ylabel("$\omega$ (rad/s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');

subplot(3,2,6);
grid on;hold on;box on;
plot(time_vec, UU{1}, 'LineWidth',2);
plot(time_vec, UU{2}, 'LineWidth',2, 'LineStyle','--', 'Color', dark_red);
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
set(gca,'TickLabelInterpreter','latex');
xlabel("Time (s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
ylabel("Torque (N$\,$m)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
h(1) = plot(NaN,NaN,'LineStyle',':', 'Color', grey, 'LineWidth',3);
h(2) = plot(NaN,NaN,'Marker','o','LineWidth',3,'LineStyle','none','Color',gold);
h(3) = plot(NaN,NaN,'LineWidth',3, 'Color',blue);
h(4) = plot(NaN,NaN,'LineWidth',2, 'LineStyle','--', 'Color', dark_red);
legend(h, "Training Trajectories", "System's equilibria", "Local at $x_0$","Lifted Predictor",'Interpreter','latex', 'Position',[0.144314696965392 0.830588763650145 0.306000308761985 0.165822227986653]);






