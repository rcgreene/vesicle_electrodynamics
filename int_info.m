function int_i = int_info(chnkr, smp)
% Compute necessary info to quickly find indefinite integrals
res = chnkr.k;
npt = chnkr.npt;
pan_num = npt/res;
A = CoefIntMat(res);
int_i = struct;

% Determine arcdensity, surface tangents, and curvature
d_x = reshape(chnkr.d(1, 1:end, 1:end), res, pan_num);
d_y = reshape(chnkr.d(2, 1:end, 1:end), res, pan_num);
d_xx = reshape(chnkr.d2(1, 1:end, 1:end), res, pan_num);
d_yy = reshape(chnkr.d2(2, 1:end, 1:end), res, pan_num);
int_i.n = get_geometry(d_x, d_y, d_xx, d_yy);
int_i.n_1 = get_geometry(d_x, d_y, d_xx, d_yy, smp.dmat1);
int_i.n_2 = get_geometry(d_x, d_y, d_xx, d_yy, smp.dmat2);
int_i.n_3 = get_geometry(d_x, d_y, d_xx, d_yy, smp.dmat3);
int_i.l = get_geometry(d_x, d_y, d_xx, d_yy, smp.lob * smp.n.u);

%Matrices for getting integral coefficients
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
int_i.avg1 = zeros(npt, 1);
int_i.avg3 = zeros(npt, 3);

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

function g_s = get_geometry(d_x, d_y, d_xx, d_yy, dmat)
    if exist("dmat", "var")
        d_x = dmat * d_x;
        d_y = dmat * d_y;
        d_xx = dmat * d_xx;
        d_yy = dmat * d_yy;
    end
    g_s = struct;
    ds = sqrt(d_x.^2 + d_y.^2);
    
    g_s.xs = [d_x(:)'; d_y(:)']./ds(:)';
    g_s.xss = [d_xx(:)'; d_yy(:)']./(ds(:).^2)';
    tmp = sum(g_s.xs .* g_s.xss, 1);
    g_s.xss = g_s.xss - g_s.xs .* tmp;
    g_s.s = ds;
end

function alt_sum = two_sum(in_mat)
    alt_sum = zeros(size(in_mat, 2), 2);
    alt_sum(1:end, 2) = sum(in_mat, 1);
    alt_sum(1:end, 1) = alt_sum(1:end, 2) - 2*sum(in_mat(2:2:end, 1:end), 1)';
end