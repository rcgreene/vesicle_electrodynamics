function ss = resample_info(n)
% Structure containing all of the necessary info for resampling
ss = {};
ss.n = {};
ss.n_1 = {};
ss.n_2 = {};
ss.n_3 = {};
[ss.n.x, ss.n.w, ss.n.u, ss.n.v] = lege.exps(n);
[ss.n_1.x, ss.n_1.w, ss.n_1.u, ss.n_1.v] = lege.exps(n-1);
[ss.n_2.x, ss.n_2.w, ss.n_2.u, ss.n_2.v] = lege.exps(n-2);
[ss.n_3.x, ss.n_3.w, ss.n_3.u, ss.n_3.v] = lege.exps(n-3);
ss.dmat = ss.n_1.v * eye(n - 1, n) * ss.n.u;
ss.umat = ss.n.v * eye(n, n-1) * ss.n.u_1;

lob_v = zeros(n, n);
lob_v(1:end, 1) = 1;
x = lobatto_nodes(n);
for i = 2:n
    lob_v(1:end, i) = lege.pol(x, i - 1);
end
ss.loba = lob_v; %resample at lobatto nodes from Legendre coefs.
end