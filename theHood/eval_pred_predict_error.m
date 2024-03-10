function errors = eval_pred_predict_error(sys_struct, X, U, cost_funcs, t)
%eval_pred_predict_error Evaluate the prediction error of the predictor

% Init
errors = 0*cost_funcs(X{1}(:,1), X{1}(:,end)); % Should be zeros of size: number of errors

% Unpack
sys = sys_struct.sys;
lift_f = sys_struct.lift_func;

Ntraj = size(X,1);
t_cut = t(1:size(U{1},2));

% For every trajectory
for ii = 1:Ntraj
    Xcurr_eval = X{ii};
    Ucurr = U{ii};
    X0 = Xcurr_eval(:,1);
    
    % Use predictor to predict the trajectory
    x_pred = lsim(sys,Ucurr,t_cut,lift_f(X0));
    x_true = Xcurr_eval';
    % figure;
    % hold on;
    % plot(Xcurr_eval');
    % plot(x_pred, 'LineStyle','--');
    % legend;

    % Vectorized: For every cost functions, compute error and update the
    % value
    errors = errors + cost_funcs(x_pred, x_true);
end

end