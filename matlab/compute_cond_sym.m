function kappa = compute_cond_sym(A_sparse)
    % COMPUTE_COND_SYM Calculates the condition number of a symmetric matrix.
    % It converts the sparse matrix to a dense format and uses the standard 
    % eigenvalue solver to reliably capture extremely small eigenvalues.
    
    % Convert to a full dense matrix
    A_dense = full(A_sparse);
    
    % Compute all eigenvalues directly
    lambdas = eig(A_dense);
    
    % Calculate the condition number as the ratio of the largest 
    % to the smallest magnitude eigenvalue
    kappa = max(abs(lambdas)) / min(abs(lambdas));
end