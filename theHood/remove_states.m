function [A_lin_reduc,B_lin_reduc] = remove_X_psi(A_lin, B_lin)

A_lin_reduc = A_lin([1,2,3,5], [1,2,3,5]);
B_lin_reduc = B_lin([1,2,3,5], :);

end
