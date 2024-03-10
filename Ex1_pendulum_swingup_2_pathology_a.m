addpath("theHood/");
addpath("theHood/ode_solvers/")

clear;
rng(10);

% IMPORTANT NOTE:
% Results depend on Matlab version (tested on MATLAB 2023b) since some
% lifting functions yield poor results and the LQR design might not be achievable 
% Also, results is heavily influenced by the sequence of random numbers.
% Nevertheless, the main conclusion is that high-order liftings are not
% always suitable.

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
Ntraj_val = 1000;   % Number of validation trajectories
AMP = 0.01;
X0_amp = 0.3*pi;
tspan = 0:Ts:Tsim;

% Final simulation
x0 = [pi;0];
time_vec = 0:Ts:1;

%% Collect Idf data
u_amp = AMP*rand(1, Ntraj);
f_u = @(t, x, traj)  0.001*randn(1) + u_amp(traj)*(- x(1,:));

X0 = [pi-X0_amp*(randn(1, Ntraj)); 1*randn(1, Ntraj)];
[X,Y,U] = collect_trajs(f_pend, f_u, X0, Ntraj, tspan, 0, "ode4");
[Xsep_idf, ~, Usep_idf] = separate_trajs(X,Y,U, Ntraj, tspan(end), Ts);

%% Identification
idx = 1;

% 1) Linearized system
A_lin = [0, 1; g/l, -gamma];
B_lin = [0; 1/(m*l^2)];
C_lin = eye(2);
systems{idx}.sys = c2d(ss(A_lin,B_lin,C_lin,0, 0), Ts);
systems{idx}.lift_func = @(X) X;
idx = idx+1;

% 2) Just sinus
lift_f1 = @(X) [X; sin(X(1,:))];
[A,B,C]= SystemID_via_EDMD(X,Y,U, lift_f1);
systems{idx}.sys = ss(A,B,C,0, Ts);
systems{idx}.lift_func = lift_f1;
idx = idx+1;

% ======================= Higher-order LLP =======================
N_lift = 100;
% ========== Splines ===========
N_splines = N_lift;
y0 = 10*(rand(1,N_splines)-0.5);
Splines_f = cell(N_splines,1);
for ii = 1:N_splines
    Splines_f{ii} = @(y) norm_vec(y(1,:) - y0(:,ii)).^2.*log(norm_vec(y(1,:) - y0(:,ii)));
end
temp = "y";
for jj = 1:N_splines
    temp = temp + " ; " + "Splines_f{"+num2str(jj)+"}(y)";
end
Lift_function = eval("@(y) [" + temp + "]");
[A,B,C]= SystemID_via_EDMD(X,Y,U, Lift_function);
systems{idx}.sys = ss(A,B,C,0, Ts);
systems{idx}.lift_func = Lift_function;
idx = idx+1;

% ========== Monomials ==========
N_polynoms = N_lift;
Polynoms_f = cell(N_polynoms,1);
for ii = 1:N_polynoms
    Polynoms_f{ii} = @(y) y(1,:).^(ii+1)./factorial(ii+1);
end
temp = "y; y(1,:)*0+1";
for jj = 1:N_polynoms
    temp = temp + " ; " + "Polynoms_f{"+num2str(jj)+"}(y)";
end
Lift_functions = eval("@(y) [" + temp + "]");
[A,B,C]= SystemID_via_EDMD(X,Y,U, Lift_functions);
systems{idx}.sys = ss(A,B,C,0, Ts);
systems{idx}.lift_func = Lift_functions;
idx = idx+1;

% ======================= RBF basis =======================
N_RBF = N_lift;
Lift_functions = cell(N_RBF, 1);

mu = 2*rand(1, N_RBF);
s = 2*rand(1, N_RBF);
RBF_f = cell(N_RBF,1);
for ii = 1:N_RBF
    RBF_f{ii} = @(y) exp( -(y(1,:) -mu(ii)).^2 ./ (2*s(ii)^2));
end
for ii = N_RBF
    temp = "y";
    for jj = 1:ii
        temp = temp + " ; " + "RBF_f{"+num2str(jj)+"}(y)";
    end
    Lift_functions{ii,1} = eval("@(y) [" + temp + "]");
    [A,B,C]= SystemID_via_EDMD(X,Y,U, Lift_functions{ii,1});
    systems{idx}.sys = ss(A,B,C,0, Ts);
    systems{idx}.lift_func = Lift_functions{ii};
    idx = idx+1;
end

%% Gather trajectories for evaluation
X0_eval = [pi-X0_amp*(randn(1, Ntraj_val)); 1*randn(1, Ntraj_val)];
u_amp = AMP*rand(1, Ntraj_val);
f_u = @(t, x, traj) 0.001*randn(1) + u_amp(traj)*(- x(1,:));

[X_eval,Y_eval,U_eval] = collect_trajs(f_pend, f_u, X0_eval, Ntraj_val, tspan, 0, "ode4");
[Xsep, Ysep, Usep] = separate_trajs(X_eval,Y_eval,U_eval, Ntraj_val, Tsim, Ts);

%% Error computation
Pred_error = cell(1, numel(systems));
Step_errors = cell(1, numel(systems));

for idx_sys = 1:numel(systems)
    % ------ Prediction Errors ------
    cost_funcs_pred = @(x, x_eval) mean(norm_vec((x - x_eval)').^2)/Ntraj_val;
    Pred_error{idx_sys} = eval_pred_predict_error(systems{idx_sys}, Xsep, Usep, cost_funcs_pred, tspan);
    
    % ------ Step Errors ------
    cost_funcs_step = @(x, u, x_n, A, B, C, lift_f) [
        mean( norm_vec(x_n - C*(A*lift_f(x) + B*u)).^2);
        mean( norm_vec(lift_f(x_n) - (A*lift_f(x) + B*u)).^2);
        ]/Ntraj_val;
    Step_errors{idx_sys} = eval_pred_step_error(X_eval, Y_eval, U_eval, cost_funcs_step, systems{idx_sys});
end

Step_errors_mat = (cell2mat(Step_errors));
Pred_error_mat = (cell2mat(Pred_error));
%% Closed-loop for every Predictor
XX = cell(1, numel(systems));
UU = cell(1, numel(systems));

% For each predictor
for predictor_idx = 1:numel(systems)
    curr_predictor = systems{predictor_idx};
    Alift = curr_predictor.sys.A;
    Blift = curr_predictor.sys.B;
    Clift = curr_predictor.sys.C;
    lift_func = curr_predictor.lift_func;
    
    % Design LQR:
    Q = blkdiag(diag([1, 1]), zeros(size(Alift,1)-2));
    R = 1000;
    [K,~,~] = dlqr(Alift,Blift,Q,R,0);
    
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

%% Compare prediction abilities
f_u = @(t, x, traj) -0.01*(x(1,:));
t_compare = 0:Ts:1;
t_cut = t_compare(1:end-1);
[X_compare,~,U_compare] = collect_trajs(f_pend, f_u, [pi;0], 1, t_compare, 0, "ode4");

XX_pred = cell(1, numel(systems));
for predictor_idx = 1:numel(systems)
    sys = systems{predictor_idx}.sys;
    lift_f = systems{predictor_idx}.lift_func;
    x_pred = lsim(sys,U_compare,t_cut,lift_f(X_compare(:,1)));
    XX_pred{predictor_idx} = x_pred';
end

all_err = [Step_errors_mat; Pred_error_mat]';

%% Generate figures
gca_fontsize = 12;
blue = [0 0.447058823529412 0.741176470588235];
dark_red = [0.635294117647059 0.0784313725490196 0.184313725490196];
gold = [0.929411764705882 0.694117647058824 0.125490196078431];
grey = [0.501960784313725 0.501960784313725 0.501960784313725];
set(0,'defaulttextinterpreter','latex');
FigHandle = figure('Position', [100, 100, 600, 600]); % nice proportions for 2-column article

subplot(3,2,1:2);
bar(all_err');
set(gca,'YScale','log')
hold on;
grid on;
box on;
ylim([10^-10, 10^6]);
Leg = legend("Local at $x_0$", "Sine", "TPS", "Polynoms", "Gaussian", 'Interpreter','latex', ...
    'Position',[0.0917219332837482 0.934246618770358 0.848633944283904 0.0472888903299966],'Orientation','horizontal');
set(gca, 'XTickLabel', {'$\epsilon_\mathrm{projected}$','$\epsilon_\mathrm{lifted}$','$e_\mathrm{prediction}$'});
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
ylabel("Error (-)")
ax=gca;
ax.XAxis.FontSize = 17;

set(gca,'TickLabelInterpreter','latex');
subplot(3,2,3);
grid on;hold on;box on;
plot(t_cut,XX_pred{1}(1,:), 'LineWidth',2);
plot(t_cut,XX_pred{2}(1,:), 'LineWidth',2);
plot(t_cut,XX_pred{3}(1,:), 'LineWidth',2);
plot(t_cut,XX_pred{4}(1,:), 'LineWidth',2);
plot(t_cut,XX_pred{5}(1,:), 'LineWidth',2);
plot(t_cut,X_compare(1,:), 'LineStyle','--', 'LineWidth',2, 'Color', 'black');
ylim([-2 4]);
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
set(gca,'TickLabelInterpreter','latex');
ylabel("$\varphi$ (rad)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
xlabel("Time (s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');

subplot(3,2,4);
grid on;hold on;box on;
plot(t_cut,XX_pred{1}(2,:), 'LineWidth',2);
plot(t_cut,XX_pred{2}(2,:), 'LineWidth',2);
plot(t_cut,XX_pred{3}(2,:), 'LineWidth',2);
plot(t_cut,XX_pred{4}(2,:), 'LineWidth',2);
plot(t_cut,XX_pred{5}(2,:), 'LineWidth',2);
plot(t_cut,X_compare(2,:), 'LineStyle','--', 'LineWidth',2, 'Color', 'black');
ylim([-20 20]);
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
set(gca,'TickLabelInterpreter','latex');
ylabel("$\omega$ (rad/s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
xlabel("Time (s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
h(1) = plot(NaN,NaN,'LineStyle','--', 'Color', 'black', 'LineWidth',2);
legend(h, "True Trajectory",  'Interpreter','LaTex', ...
    'Position',[0.698448136005573 0.607837991750711 0.251866869721804 0.0354666677474973]);

subplot(3,2,5:6);
grid on;hold on;box on;
% title("Closed-loop performance")
plot(time_vec, XX{1}(1,:), 'LineWidth',2);
plot(time_vec, XX{2}(1,:), 'LineWidth',2, 'LineStyle','--');
plot(time_vec, XX{3}(1,:), 'LineWidth',2, 'LineStyle',':');
plot(time_vec, XX{4}(1,:), 'LineWidth',2, 'LineStyle','-.');
plot(time_vec, XX{5}(1,:), 'LineWidth',2, 'LineStyle','--');
ylim([-pi 1.5*pi]);
set(gca,'fontsize',gca_fontsize); % set size of the labels at xlabel and ylabel
set(gca,'TickLabelInterpreter','latex');
ylabel("$\varphi$ (rad)",'FontSize',gca_fontsize, 'Interpreter','LaTex');
xlabel("Time (s)",'FontSize',gca_fontsize, 'Interpreter','LaTex');





