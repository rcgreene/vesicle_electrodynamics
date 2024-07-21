function STau = STau_kernel(mu)

STau = kernel();
STau.name = 'stokes';
STau.params.mu = mu;

STau.type = 'stt';
STau.eval = @(s, t) k_func(mu, s, t);
STau.opdims = [2, 2];
STau.sing = 'pv';

end

function varargout = k_func(mu, s, t)

    rx = t.r(1, :).' - s.r(1, :);
    ry = t.r(2, :).' - s.r(2, :);
    r2 = rx.^2 + ry.^2;
    [Nt, Ns] = size(rx);

    r4 = r2.^2;
    nx = t.n(1, :)';
    ny = t.n(2, :)';
    tx = -rx.*ny;
    ty = ry.*nx;
    sq_diff = (ry + rx).*(ry - rx);
    cf = 1/(4*pi*mu);
    Kxx = cf * ((tx .* sq_diff - ty .* (3*rx.^2 + ry.^2)) ./r4);
    Kyy = cf * ((-ty .* sq_diff - tx .* (3*ry.^2 + rx.^2)) ./r4);
    Kxy = -cf * (((ny .* ry + nx .* rx) .* sq_diff) ./r4);

    if (nargout == 1)
        K = zeros(2*Nt, 2*Ns);
        K(1:2:end, 1:2:end) = Kxx;
        K(2:2:end, 1:2:end) = Kxy;
        K(1:2:end, 2:2:end) = Kxy;
        K(2:2:end, 2:2:end) = Kyy;
        varargout = {K};
    else
        varargout = {Kxx, Kxy, Kyy};
    end
end