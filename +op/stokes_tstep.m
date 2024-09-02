function [sys_mat, sys_dmat] = stokes_tstep(chnkr, stt_mat, smp, int_i, dt)

N = chnkr.npt;
res = chnkr.k;
p_num = N/res;

[B, I_op] = op.bend_eval(chnkr, stt_mat, smp, int_i);
sys_mat = (I_op - dt*B);
dmat_vec = zeros(2*res, 2*res);
dmat_vec(7:2:end, 1:2:end) = smp.l2.dmat(2:end, 1:end);
dmat_vec(8:2:end, 2:2:end) = smp.l2.dmat(2:end, 1:end);
sys_dmat = kron(speye(p_num), dmat_vec);
sys_mat = kron(speye(p_num), dmat_vec) * sys_mat;

sys_mat(1:3, 1:2:end) = int_i.avg3';
sys_mat(4:6, 2:2:end) = int_i.avg3';

for i = 2:p_num
    j = chnkr.adj(2, i);
    for k = 1:3
        sys_mat((i-1)*2*res + 2*k - 1, (i-1)*2*res + 1:2:i*2*res) = int_i.cont3(1:end, 2, k, i);
        sys_mat((i-1)*2*res + 2*k, (i-1)*2*res + 2:2:i*2*res) = int_i.cont3(1:end, 2, k, i);
        sys_mat((i-1)*2*res + 2*k - 1, (j-1)*2*res + 1:2:j*2*res) = -int_i.cont3(1:end, 1, k, j);
        sys_mat((i-1)*2*res + 2*k, (j-1)*2*res + 2:2:j*2*res) = -int_i.cont3(1:end, 1, k, j);
    end
end

end