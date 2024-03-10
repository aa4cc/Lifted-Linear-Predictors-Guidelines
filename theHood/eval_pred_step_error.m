function errors = eval_pred_step_error(X_eval, Y_eval, U_eval, cost_funcs, system_struct)
    
    A = system_struct.sys.A;
    B = system_struct.sys.B;
    C = system_struct.sys.C;
    lift_f = system_struct.lift_func;

    errors = cost_funcs(X_eval, U_eval, Y_eval, A, B, C, lift_f);
end