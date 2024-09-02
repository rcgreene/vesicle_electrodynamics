function D = D_bcs(chnkr, int_i, smp, dv_mat)

    N = chnkr.npt;
    p = chnkr.npt/chnkr.k;
    i = 1:2*N;
    j = floor((i + 1)./2);
    xs_n_1 = zeros(2, chnkr.k, p);
    %xs_n_1(1:2, 2:end, 1:end) = reshape(int_i.n_1.xs, 2, chnkr.k - 1, p);
    xs_n_1(1:2, 1:end, 1:end) = reshape(int_i.l.xs, 2, chnkr.k, p);
    xs_n_1(1:2, 1, 1:end) = 0;
    s_n_1 = zeros(chnkr.k, p);
    s_n_1(1:end, 1:end) = 1./int_i.l.s;
    s_n_1(1, 1:end) = 0;
    xs_mat = sparse(j, i, xs_n_1(:));
%    s_n_1 = zeros(chnkr.k, p);
%    s_n_1(2:end, 1:end) = 1./chnkr.n_1.s;
    tmp = zeros(2*chnkr.k, 2*chnkr.k);
    tmp(1:2:end, 1:2:end) = smp.lob * smp.n.u;
    tmp(2:2:end, 2:2:end) = smp.lob * smp.n.u;
    tmp(1:2, 1:end, 1:end) = 0;
    %tmp = zeros(chnkr.k, chnkr.k);
    %tmp(2:end, 1:end) = smp.n_1.v * eye(chnkr.k - 1, chnkr.k) * smp.n.u;

%    D = kron(speye(p), tmp) * sparse(1:N, 1:N, 1./int_i.n.s(:)) ...
%        * xs_mat * kron(speye(p), dv_mat);

    D = sparse(1:N, 1:N, s_n_1(:)) ...
        * xs_mat * kron(speye(p), tmp)  * kron(speye(p), dv_mat);

end