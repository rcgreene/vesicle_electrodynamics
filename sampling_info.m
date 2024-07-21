function smp = sampling_info(n)
% Structure containing all of the necessary info for resampling
smp = {};
smp.n = {};
smp.n_1 = {};
smp.n_2 = {};
smp.n_3 = {};
[smp.n.x, smp.n.w, smp.n.u, smp.n.v] = lege.exps(n);
[smp.n_1.x, smp.n_1.w, smp.n_1.u, smp.n_1.v] = lege.exps(n-1);
[smp.n_2.x, smp.n_2.w, smp.n_2.u, smp.n_2.v] = lege.exps(n-2);
[smp.n_3.x, sMp.n_3.w, smp.n_3.u, smp.n_3.v] = lege.exps(n-3);
smp.dmat1 = smp.n_1.v * eye(n - 1, n) * smp.n.u;
smp.umat1 = smp.n.v * eye(n, n-1) * smp.n_1.u;
smp.umat3 = smp.n.v * eye(n, n-3) * smp.n_3.u;

lob_v = zeros(n, n);
lob_v(1:end, 1) = 1;
t = sin(linspace(-pi/2, pi/2, 2*n));
x = t(1:2:end);
%x = lobatto_nodes(n);
for i = 2:n
    lob_v(1:end, i) = lege.pol(x, i - 1);
end
lob_u = inv(lob_v);
smp.lob_u = lob_u;
smp.lob = lob_v; %resample at lobatto nodes from Legendre coefs.
lob_v2 = zeros(n - 2, n);
lob_v2(1:end, 1) = 1;
x = lobatto_nodes(n - 2);
for i = 2:n
    lob_v2(1:end, i) = lege.pol(x, i - 1);
end
lob_u2 = pinv(lob_v2);
smp.lob_u2 = lob_u2;
smp.lob_v2 = lob_v2(1:n-2, 1:n-2);
smp.l2.dmat = lob_v2 * smp.n.u;
% Set up the down sampling matrices for smaller matrices
alias_3 = AliasingMat(n);
tmp_mat = eye(n - 3, n);
tmp_mat(end, end - 1) = alias_3(1,1);
tmp_mat(end - 1, end) = alias_3(2,2);
smp.dmat3 = smp.n_3.v*tmp_mat*smp.n.u;
tmp_mat = eye(n - 2, n);
tmp_mat(end, end) = -(n - 2)/(n - 1);
smp.dmat2 = smp.n_2.v*tmp_mat*smp.n.u;

end

function output = AliasingMat(N)
    % Determine the appropriate aliasing relation for legendre polynomials
    % of order N - p + x onto legendre polynomials of order N - p - x with
    % an N - p point grid. (This uses a standard recurrence relation)
    % p is only ever 3 or 1, which trivializes the function, so I've
    % dummied it
    p = 3;
    output = zeros(p - 1, p - 1);
    tmp_val = -(N - p)/(N - p + 1);
    output(1, 1) = tmp_val;
    tmp_val = tmp_val *(2*(N - p) + 3)/(N - p + 2);
    tmp_val = tmp_val * (N - p - 1)/(2*(N - p) - 1);
    output(2, 2) = tmp_val;
end