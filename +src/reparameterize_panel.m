function [r_slice_out] = reparameterize_panel(r_slice_in)

    [~, n] = size(r_slice_in);
    r_slice_out = r_slice_in;
    [leg_roots, ~, u, v] = lege.exps(n);
    tmp = zeros(n, n);
    tmp(1:end - 1, 1:end) = lege.derpol(eye(n));
    du = v * tmp * u;
    s_func = @(r) sqrt( sum((r * du').^2, 1));

    tmp = spdiags([-1./(1:2:(2*n - 1))' 1./(1:2:(2*n - 1))'], [1 -1], n, n);
    int_u = v * tmp * u;

    err = 1;
    while(err > 1e-5)
        s = s_func(r_slice_out)';
        arc_fac = (u(1, 1:end) * s);
        s = s ./ arc_fac;
        s_int = int_u * s;
        resid = leg_roots - s_int;
        t = leg_roots + (leg_roots - s_int)./s;
        coefs = r_slice_out * u';
        r_slice_out = [lege.exev(t, coefs(1, 1:end)')'; lege.exev(t, coefs(2, 1:end)')'];
        err = norm(resid);
    end

end