function mu = mu_update(int_i, mu, smp)

    mu_r = reshape(mu, 2, chnkr.k, chnkr.npt/chnkr.k);

    x_ttt = mu_r(1, 4:end, 1:end);
    y_ttt = mu_r(2, 4:end, 1:end);

    % This would be so much easier if einstein summation existed in matlab

    x_tt = int_i.n_3.xss(1,:);
    y_tt = int_i.n_3.xss(2, :);
    x_tt = x_tt(:);
    y_tt = y_tt(:);
    x_s = int_i.n_3.xs(1, 1:end)';
    y_s = int_i.n_3.xs(2, 1:end)';
    x_ss = int_i.n_3.xss(1, 1:end)';
    y_ss = int_i.n_3.xss(2, 1:end)';
    s_p = x_tt.*x_s + y_tt .* y_s;


    x_p = x_ttt(:)./int_i.n_3.s(:).^3 - 2*s_p(:).*x_ss(:)./int_i.n_3.s(:).^2;
    y_p = y_ttt(:)./int_i.n_3.s(:).^3 - 2*s_p(:).*y_ss(:)./int_i.n_3.s(:).^2;
    s_pp = x_p .* x_s + y_p .* y_s;
    x_p = x_p - s_pp .* x_s;
    y_p = y_p - s_pp .* y_s;
    
    s_pp = -(x_tt .* x_s + y_tt .* y_s)./(int_i.n_3.s(:).^3);
    x_p = x_p + x_ss.*s_pp;
    y_p = y_p + y_ss.*s_pp;
    s_pp = -(x_tt .* x_ss + y_tt .* y_ss)./(int_i.n_3.s(:).^3);
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
end