bell_opts = [];
bell_opts.type = 'prolate';
%addpath("/home/connor/projects/funextend/funextend2d/matlab/")
%addpath("/home/connor/bin/finufft/matlab/")
cparams = [];
cparams.eps = 1.0e-12;
cparams.nover = 0;
cparams.maxchunklen = .5;

K_type = "combined";
do_reparam = false;

panel_res = 16; 
pref = []; 
pref.k = panel_res;


k_ss = src.STS_kernel(1);
k_st = src.STau_kernel(1);
k_s  = kernel("stokes", "s", 1);

    arclen = dictionary();
    volume = dictionary();

T_end = 10.;
dt_list = [];
arclen_err = [];
volume_err = [];

for iter = 100:100

chnkr_in = chunkerfunc(@(t) [(1 - .4*cos(2*t))'.*(cos(t))'; 2*(1 - .4*cos(2*t))'.*(sin(t))'],cparams,pref); 
N = chnkr_in.npt;
pan_num = chnkr_in.npt/panel_res;
if do_reparam
    for i = 1:pan_num
        chnkr_in.r(1:2, 1:end, i) = src.reparameterize_panel(chnkr_in.r(1:2, 1:end, i));
    end
    chnkr_in = update_chunker(chnkr_in)
end
%Get relevent legendre info
smp = sampling_info(panel_res);
dv = smp.n.v * eye(panel_res, panel_res-1) * lege.derpol(smp.n.u);
dv_mat = zeros(2*panel_res, 2*panel_res);
dv_mat(1:2:end, 1:2:end) = dv;
dv_mat(2:2:end, 2:2:end) = dv;

mu = get_mu(chnkr_in, smp);
int_i = int_info(chnkr_in, smp);
int_j = int_info(chnkr_in, smp);
plot(chnkr_in)
hold on
[v_1, v_2] = conserved_qs(chnkr_in);

dt = T_end/(20*iter)

dt_list = [dt_list dt];

arclen{dt} = v_1;
volume{dt} = v_2;



for i = 1:(20*iter)

    S_mat  = chunkermat(chnkr_in, k_s);
    if K_type == "combined"
        ST_mat = chunkermat(chnkr_in, k_st);
        SS_mat = chunkermat(chnkr_in, k_ss);
    else
        [D_1, D_2] = op.D_no_bcs(chnkr_in, int_i, dv_mat);
        ST_mat = D_2 * S_mat;
        SS_mat = S_mat * D_2;

    end
    tic;
    [tens_mat, shift_mat] = op.tens_solve(chnkr_in, ST_mat, smp, int_i);
    toc
    [sys_mat, sys_dmat]   = op.stokes_tstep(chnkr_in, SS_mat, smp, int_i, dt);
    K = op.bend_eval(chnkr_in, SS_mat, smp, int_i);
    tic;
    T = op.tens_eval(chnkr_in, S_mat, smp, int_i);
    toc
    tic;
    D = op.D_bcs(chnkr_in, int_i, smp, dv_mat);
    toc
    rhs_1 = -D*K*mu;
    tic;
    sigma = gmres(tens_mat, rhs_1, [], 1e-12, 1000);
    %norm(D*(K*mu + T*sigma))
    %norm(K*mu + T*sigma)
    toc
    rhs_2 = dt*(T*sigma) + chnkr_in.r(:);
    %rhs_2(:) = ST_mat * rhs_2;
    rhs_2 = sys_dmat * rhs_2;
    tic;
    mu_new = gmres(sys_mat, rhs_2, [], 1e-12, 3000, [], [], mu);
    toc
    %[a, b, c] = update_int_info(int_i, chnkr_in, mu_new, smp);
    [int_i, r, mu] = update_int_info(int_i, chnkr_in, mu_new, smp);
    norm(chnkr_in.r(:) - r(:))
    norm(D * (chnkr_in.r(:) - r(:)))
    chnkr_in.r(:) = r(:);
    chnkr_in = update_chunker(chnkr_in, chnkr_in.r);
    [v_1, v_2] = conserved_qs(chnkr_in);
    arclen{dt} = [arclen{dt} v_1];
    volume{dt} = [volume{dt} v_2];
    if (mod(i, 1) == 0)
        plot(chnkr_in)
        hold on
        pause(.01)
    end
end

arclen_err = [arclen_err, arclen{dt}(1) - arclen{dt}(end)];
volume_err = [volume_err, volume{dt}(1) - volume{dt}(end)];

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