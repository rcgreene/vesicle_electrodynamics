function [D_1, D_2] = D_no_bcs(chnkr, int_i, dv_mat)

    N = chnkr.npt;
    p = chnkr.npt/chnkr.k;
    i = 1:2*N;
    j = floor((i + 1)./2);
    %xs_n_1(1:2, 2:end, 1:end) = reshape(int_i.n_1.xs, 2, chnkr.k - 1, p);
    xs_n = reshape(int_i.n.xs, 2, chnkr.k, p);
    s_n = 1./int_i.n.s;
    xs_mat = sparse(j, i, xs_n(:));

    D_1 = sparse(1:N, 1:N, s_n(:)) ...
        * xs_mat  * kron(speye(p), dv_mat);

    D_2 = (sparse(1:2:(2*N), 1:2:(2*N), s_n(:), 2*N, 2*N) + ...
           sparse(2:2:(2*N), 2:2:(2*N), s_n(:), 2*N, 2*N)) * kron(speye(p), dv_mat);
end