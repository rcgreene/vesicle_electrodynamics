bell_opts = [];
bell_opts.type = 'prolate';
%addpath("/home/connor/projects/funextend/funextend2d/matlab/")
%addpath("/home/connor/bin/finufft/matlab/")
cparams = [];
cparams.eps = 1.0e-11;
cparams.nover = 0;
cparams.maxchunklen = .1;

panel_res = 16; 
pref = []; 
pref.k = panel_res;

chnkr_in = chunkerfunc(@(t) [(1 - .29*cos(2*t) - .15*cos(6*t))'.*(cos(t))'; 2*(1 - .29*cos(2*t) - .15*cos(6*t))'.*(sin(t))'],cparams,pref); 
N = chnkr_in.npt;
pan_num = chnkr_in.npt/panel_res;
%Get relevent legendre info
smp = sampling_info(panel_res);
dv = smp.n.v * eye(panel_res, panel_res-1) * lege.derpol(smp.n.u);
dv_mat = zeros(2*panel_res, 2*panel_res);
dv_mat(1:2:end, 1:2:end) = dv;
dv_mat(2:2:end, 2:2:end) = dv;

mu = get_mu(chnkr_in, smp);
int_i = int_info(chnkr_in, smp);
int_j = int_info(chnkr_in, smp);
k_ss = STS_kernel(1);
k_st = STau_kernel(1);
k_s  = kernel("stokes", "s", 1);
plot(chnkr_in)
hold on
dt = .02;

for i = 1:20

    S_mat = chunkermat(chnkr_in, k_s);
    ST_mat = chunkermat(chnkr_in, k_st);
    SS_mat = chunkermat(chnkr_in, k_ss);
    tic;
    [tens_mat, shift_mat] = make_tension_matrix(chnkr_in, ST_mat, smp, int_i);
    toc
    [sys_mat, sys_dmat]= stokes_system_matrix(chnkr_in, SS_mat, smp, int_i, dt);
    K = make_bending_matrix(chnkr_in, SS_mat, smp, int_i);
    tic;
    T = tension_eval_matrix(chnkr_in, S_mat, smp, int_i);
    toc
    tic;
    D = make_D_matrix(chnkr_in, int_i, smp, dv_mat);
    toc
    rhs_1 = -shift_mat*D*K*mu;
    tic;
    sigma = gmres(tens_mat, rhs_1, [], 1e-6, 1000);
    norm(D*(K*mu + T*sigma))
    norm(K*mu + T*sigma)
    toc
    rhs_2 = dt*(T*sigma + K*mu) + chnkr_in.r(:);
    %rhs_2(:) = ST_mat * rhs_2;
    rhs_2 = sys_dmat * rhs_2;
    tic;
    mu_new = gmres(sys_mat, rhs_2, [], 1e-6, 3000, [], [], mu);
    toc
    %[a, b, c] = update_int_info(int_i, chnkr_in, mu_new, smp);
    [int_i, r, mu] = update_int_info(int_i, chnkr_in, mu_new, smp);
    chnkr_in.r(:) = r(:);
    chnkr_in = update_chunker(chnkr_in, chnkr_in.r);
    if (mod(i, 1) == 0)
        plot(chnkr_in)
        hold on
        pause(.01)
    end
end

%r = mu_new;
%for i = 1:pan_num
%    j = (i - 1)*panel_res*2 + 1;
%    k = i*panel_res*2;
%    r(j + 4:2:k) = smp.n_2.v * int_i.int3(1:end, 1:end, i) * r(j+4:2:k);
%    r(j + 5:2:k) = smp.n_2.v * int_i.int3(1:end, 1:end, i) * r(j+5:2:k);
%    r(j + 2:2:k) = smp.n_1.v * int_i.int2(1:end, 1:end, i) * r(j+2:2:k);
%    r(j + 3:2:k) = smp.n_1.v * int_i.int2(1:end, 1:end, i) * r(j+3:2:k);
%    r(j:2:k) = smp.n.v * int_i.int1(1:end, 1:end, i) * r(j:2:k);
%    r(j + 1:2:k) = smp.n.v * int_i.int1(1:end, 1:end, i) * r(j+1:2:k);
%end

%chnkr_in.r(:) = r;