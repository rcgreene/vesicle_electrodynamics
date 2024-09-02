function [bend_mat, int_mat_2] = bend_eval(chnkr, stt_kernel, smp, int_i)

N = chnkr.npt;
res = chnkr.k;
p_num = N/res;

a = 0:(p_num*res^2 - 1);
i = (res)*floor(a/(res^2)) + mod(a, res) + 1;
j = (res)*floor(a/(res^2)) + floor(mod(a, res.^2)/res) + 1;

int3 = zeros(res, res, p_num);
int2 = zeros(res, res, p_num);
int3(1, 1, 1:end) = 1;
int3(2, 2, 1:end) = 1;
int2(1, 1, 1:end) = 1;
int3(3:end, 3:end, 1:end) = int_i.int3;
int2(2:end, 2:end, 1:end) = int_i.int2;
v_2 = eye(res);
v_2(3:end, 3:end) = smp.n_2.v;
v_1 = eye(res);
v_1(2:end, 2:end) = smp.n_1.v;

% integrate three times along arc
int_mat = kron(eye(p_num), smp.n.v) * sparse(i, j, int_i.int1(:)) * ...
    kron(eye(p_num), v_1) * sparse(i, j, int2(:)) * ...
    kron(eye(p_num), v_2) * sparse(i, j, int3(:));
int_mat_2 = zeros(2*N, 2*N);
int_mat_2(1:2:end, 1:2:end) = int_mat;
int_mat_2(2:2:end, 2:2:end) = int_mat;
% apply kernel
tmp = zeros(2*res, 2*res);
tmp(1:2:end, 7:2:end) = smp.umat3;
tmp(2:2:end, 8:2:end) = smp.umat3;

bend_mat = -stt_kernel * kron(eye(p_num), tmp);
end