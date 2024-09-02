function output = CoefIntMat(N)
    % Maps N-p coefficients to the pth order indefinite integral with 0 
    % constants of integration.
    tmp_vec = 1./(2.*(0:(N)) + 1);
    tmp_mat = spdiags([-(tmp_vec(1:end-2))' tmp_vec(1:end-2)'], [2 0], N-1, N-1);
    output = tmp_mat;
end
