function outmat = tension_eval_matrix(chnkr, s_kernel, smp, int_i)

    N = chnkr.npt;
    res = chnkr.k;
    p_num = N/res;
    i = 1:2*N;
    j = floor((i + 1)./2);

    xs_mat = sparse(i, j, int_i.n.xs(:));
    xss_mat = sparse(i, j, int_i.n.xss(:));

    intmat = zeros(N, N);
    umat = zeros(res, res);
    umat(1:end, 2:end) = smp.umat1;
    upmat = kron(speye(p_num), umat);
    for i = 1:p_num
        intmat((i - 1)*res + 1:i*res, (i-1)*res+1:i*res) ...
            = smp.n.v  * int_i.int1(1:end,1:end,i);
    end
    outmat = s_kernel * (xss_mat * intmat + xs_mat * upmat);
end