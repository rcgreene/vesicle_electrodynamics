function [int_i, r, mu] = update_int_info(int_i, chnkr, mu, smp)
    n = chnkr.k;

    mu_r = reshape(mu, 2, chnkr.k, chnkr.npt/chnkr.k);
    
    x_tt = zeros(n-2, chnkr.npt/n);
    y_tt = zeros(n-2, chnkr.npt/n);
    x_t  = zeros(n-1, chnkr.npt/n);
    y_t  = zeros(n-1, chnkr.npt/n);
    r    = zeros(2, n, chnkr.npt/n);
    for i = 1:(chnkr.npt/chnkr.k)
        x_tt(1:end, i) = smp.n_2.v * squeeze(int_i.int3(1:end, 1:end, i)) * squeeze(mu_r(1, 3:end, i))';
        y_tt(1:end, i) = smp.n_2.v * squeeze(int_i.int3(1:end, 1:end, i)) * squeeze(mu_r(2, 3:end, i))';
        mu_r(1, 3:end, i) = x_tt(1:end, i);
        mu_r(2, 3:end, i) = y_tt(1:end, i);
        x_t(1:end, i) = smp.n_1.v * squeeze(int_i.int2(1:end, 1:end, i)) * squeeze(mu_r(1, 2:end, i))';
        y_t(1:end, i) = smp.n_1.v * squeeze(int_i.int2(1:end, 1:end, i)) * squeeze(mu_r(2, 2:end, i))';
        mu_r(1, 2:end, i) = x_t(1:end, i);
        mu_r(2, 2:end, i) = y_t(1:end, i);
    
        r(1, 1:end, i) = smp.n.v * squeeze(int_i.int1(1:end, 1:end, i)) * squeeze(mu_r(1, 1:end, i))';
        r(2, 1:end, i) = smp.n.v * squeeze(int_i.int1(1:end, 1:end, i)) * squeeze(mu_r(2, 1:end, i))';
    end


%    test_mat = 1 - x_s.^2 - y_s.^2;
%    max(test_mat(:))

    t_1 = smp.n.v * eye(n, n-1) * smp.n_1.u;
    t_2 = smp.n.v * eye(n, n-2) * smp.n_2.u;
    int_i.n = get_geometry(x_t, y_t, x_tt, y_tt, int_i.n, t_1, t_2);
    t_1 = smp.lob * eye(n, n-1)*smp.n_1.u;
    t_2 = smp.lob * eye(n, n-2)*smp.n_2.u;
    int_i.l = get_geometry(x_t, y_t, x_tt, y_tt, int_i.l, t_1, t_2);
    t_2 = smp.n_1.v * eye(n-1, n-2)*smp.n_2.u;
    int_i.n_1 = get_geometry(x_t, y_t, x_tt, y_tt, int_i.n_1, eye(n - 1), t_2);
    t_1 = smp.n_2.v * eye(n-2, n-1) * smp.n_1.u;
    int_i.n_2 = get_geometry(x_t, y_t, x_tt, y_tt, int_i.n_2, t_1, eye(n-2));
    t_1 = smp.dmat3 * smp.n.v * eye(n, n-1) * smp.n_1.u;
    t_2 = smp.n_3.v * eye(n-3, n-2) * smp.n_2.u;
    int_i.n_3 = get_geometry(x_t, y_t, x_tt, y_tt, int_i.n_3, t_1, t_2);


    mu_r = reshape(mu, 2, chnkr.k, chnkr.npt/chnkr.k);

    x_ttt = mu_r(1, 4:end, 1:end);
    y_ttt = mu_r(2, 4:end, 1:end);

    % This would be so much easier if einstein summation existed in matlab
    x_tt = t_2 * x_tt;
    y_tt = t_2 * y_tt;
    x_tt = x_tt(:);
    y_tt = y_tt(:);
    x_s = int_i.n_3.xs(1, 1:end)';
    y_s = int_i.n_3.xs(2, 1:end)';
    x_ss = int_i.n_3.xss(1, 1:end)';
    y_ss = int_i.n_3.xss(2, 1:end)';
    s_p = x_tt.*x_s + y_tt .* y_s;
    [abs(s_p), (abs(s_p) ./ sqrt(y_tt.^2 + x_tt.^2)), sqrt(x_ss.^2 + y_ss.^2), x_s.*x_ss + y_s.*y_ss];

    x_p = x_ttt(:)./int_i.n_3.ds(:).^3 - 2*s_p(:).*x_ss(:)./int_i.n_3.ds(:).^2;
    y_p = y_ttt(:)./int_i.n_3.ds(:).^3 - 2*s_p(:).*y_ss(:)./int_i.n_3.ds(:).^2;
    s_pp = x_p .* x_s + y_p .* y_s;
    x_p = x_p - s_pp .* x_s;
    y_p = y_p - s_pp .* y_s;
    
    s_pp = -(x_tt .* x_s + y_tt .* y_s)./(int_i.n_3.ds(:).^3);
    x_p = x_p + x_ss.*s_pp;
    y_p = y_p + y_ss.*s_pp;
    s_pp = -(x_tt .* x_ss + y_tt .* y_ss)./(int_i.n_3.ds(:).^3);
    x_p = x_p + x_s .* s_pp;
    y_p = y_p + y_s .* s_pp;

    mu = reshape(mu, 2, n, chnkr.npt/chnkr.k);
    mu(1, 4:end, 1:end) = reshape(x_p, 1, n - 3, chnkr.npt/chnkr.k);
    mu(2, 4:end, 1:end) = reshape(y_p, 1, n - 3, chnkr.npt/chnkr.k);

    % Extract the constant terms for mu.
    mu(1, 3, 1:end) = smp.n_2.w' * reshape(int_i.n_2.xss(1,:), n-2, chnkr.npt/n);
    mu(2, 3, 1:end) = smp.n_2.w' * reshape(int_i.n_2.xss(2,:), n-2, chnkr.npt/n);
    mu(1, 2, 1:end) = smp.n_1.w' * reshape(int_i.n_1.xs(1,:), n-1, chnkr.npt/n);
    mu(2, 2, 1:end) = smp.n_1.w' * reshape(int_i.n_1.xs(2,:), n-1, chnkr.npt/n);
    mu(1, 1, 1:end) = smp.n.w' * reshape(r(1,:), n, chnkr.npt/n);
    mu(2, 1, 1:end) = smp.n.w' * reshape(r(2,:), n, chnkr.npt/n);
    mu = mu(:);

    pan_num = chnkr.npt/chnkr.k;

% Lazy copy, should refactor later
%Matrices for getting integral coefficients
    res = chnkr.k;
    A = src.CoefIntMat(res);
    int_i.int3 = zeros(res - 2, res - 2, pan_num);
    int_i.int2 = zeros(res - 1, res - 1, pan_num);
    int_i.int1 = zeros(res, res, pan_num);
    for i = 1:pan_num
        int_i.int3(1, 1, i) = 1;
        int_i.int2(1, 1, i) = 1;
        int_i.int1(1, 1, i) = 1;
        int_i.int3(2:end, 2:end, i) = A(1:end - 2, 1:end - 2) * smp.n_3.u * diag(int_i.n_3.s(1:end, i));
        int_i.int2(2:end, 2:end, i) = A(1:end - 1, 1:end - 1) * smp.n_2.u * diag(int_i.n_2.s(1:end, i));
        int_i.int1(2:end, 2:end, i) = A * smp.n_1.u * diag(int_i.n_1.s(1:end, i));
    end

    int_i.cont1 = zeros(res, 2, pan_num);
    int_i.cont3 = zeros(res, 2, 3, pan_num);
    int_i.avg1 = zeros(chnkr.npt, 1);
    int_i.avg3 = zeros(chnkr.npt, 3);

    for i = 1:pan_num
        int_i.cont1(1:end, 1:2, i) = two_sum(int_i.int1(1:end, 1:end, i));
        int_i.avg1((i - 1)*res + 1:i*res) = (int_i.cont1(1:end, 2, i) - int_i.cont1(1:end, 1, i))/2;

        int_mat = zeros(res, res);
        int_mat(3:end, 3:end) = int_i.int3(1:end, 1:end, i);
        int_i.cont3(1:end, 1:2, 3, i) = two_sum(int_mat(3:end, 1:end));
        int_mat(3:end, 3:end) = smp.n_2.v * int_mat(3:end, 3:end);
        int_mat(2, 2) = 1;
        int_mat(2:end, 2:end) = int_i.int2(1:end, 1:end, i) * int_mat(2:end, 2:end);
        int_i.cont3(1:end, 1:2, 2, i) = two_sum(int_mat(2:end, 1:end));
        int_mat(2:end, 2:end) = smp.n_1.v * int_mat(2:end, 2:end);
        int_mat(1, 1) = 1;
        int_mat(1:end, 1:end) = int_i.int1(1:end, 1:end, i) * int_mat;
        int_i.cont3(1:end, 1:2, 1, i) = two_sum(int_mat);

        int_i.avg3((i - 1)*res + 1:i*res, 1:3) = ...
            (reshape(int_i.cont3(1:end, 2, 1:3, i) - int_i.cont3(1:end, 1, 1:3, i), res, 3))/2;
    end
end
function alt_sum = two_sum(in_mat)
    alt_sum = zeros(size(in_mat, 2), 2);
    alt_sum(1:end, 2) = sum(in_mat, 1);
    alt_sum(1:end, 1) = alt_sum(1:end, 2) - 2*sum(in_mat(2:2:end, 1:end), 1)';
end
function g_s = get_geometry(d_x, d_y, d_xx, d_yy, g_p, dmat1, dmat2)
    
    d_x = dmat1 * d_x;
    d_y = dmat1 * d_y;
    d_xx = dmat2 * d_xx;
    d_yy = dmat2 * d_yy;
    g_s = struct;
    g_s.ds = sqrt(d_x.^2 + d_y.^2);
    g_s.s = g_s.ds .* g_p.s;
    d_x = d_x(:)./g_s.ds(:);
    d_y = d_y(:)./g_s.ds(:);
    g_s.xs = [d_x(:)'; d_y(:)'];
    g_s.xss = [d_xx(:)'; d_yy(:)']./(g_s.ds(:).^2)';
    s_p = (d_x(:) .* d_xx(:) + d_y(:) .* d_yy(:))./(g_s.ds(:).^2);
    g_s.xss = g_s.xss - [s_p'.*d_x(:)'; s_p'.*d_y(:)'];
end