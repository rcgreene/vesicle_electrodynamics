bell_opts = [];
bell_opts.type = 'prolate';
%addpath("/home/connor/projects/funextend/funextend2d/matlab/")
%addpath("/home/connor/bin/finufft/matlab/")
cparams = [];
cparams.eps = 1.0e-12;
cparams.nover = 0;
cparams.maxchunklen = .1;

panel_res = 16; 
pref = []; 
pref.k = panel_res;

chnkr_in = chunkerfunc(@(t) [(1 + 0*cos(5*t))'.*(cos(t))'; (1 + 0*cos(5*t))'.*(sin(t))'],cparams,pref); 
N = chnkr_in.npt;
pan_num = chnkr_in.npt/panel_res;
%Get relevent legendre info
smp = sampling_info(panel_res);
int_i = int_info(chnkr_in, smp);

%[leg_x, leg_w, leg_u, leg_v] = lege.exps(panel_res);
%[leg_x_1, leg_w_1, leg_u_1, leg_v_1] = lege.exps(panel_res - 1);
%[leg_x_2, leg_w_2, leg_u_2, leg_v_2] = lege.exps(panel_res - 2);
%[leg_x_3, leg_w_3, leg_u_3, leg_v_3] = lege.exps(panel_res - 3);


%tmp_mat = eye(panel_res - 3, panel_res);
%tmp_mat(end, end - 1) = alias_3(1,1);
%tmp_mat(end - 1, end) = alias_3(2,2);
%down_mat_3 = leg_v_3*tmp_mat*leg_u;
%tmp_mat = eye(panel_res - 2, panel_res);
%tmp_mat(end, end) = -(panel_res - 2)/(panel_res - 1);
%down_mat_2 = leg_v_2*tmp_mat*leg_u;
%down_mat_1 = leg_v_1*eye(panel_res-1, panel_res)*leg_u;
%dmat_1_comp = leg_v * eye(panel_res, panel_res - 1)*leg_u_1;
%s = arclengthdens(chnkr_in);
%s_3 = down_mat_3 * s;
%s_2 = down_mat_2 * s;
%s_1 = down_mat_1 * s;

%int_info = struct;
%int_info.s = arclengthdens(chnkr_in);
%int_info.s_3 = smp.dmat3 * s;
%int_info.s_2 = smp.dmat2 * s;
%int_info.s_1 = smp.dmat1 * s;

%setup a differentiation matrix
%d_coef = eye(panel_res + 1);
%for i = 1:panel_res + 1
%    d_coef(1:end-1, i) = lege.derpol(d_coef(i,:)');
%end
%d_coef(end, end) = 0;
%bc_tmp = zeros(panel_res + 1, panel_res);
%bc_tmp(1:panel_res, 1:panel_res) = eye(panel_res);
%bc_tmp(panel_res + 1, 1:2:end) = -1;
%d_leg = smp.n.v * d_coef(1:panel_res, 1:panel_res) * smp.n.u;

%quickly derive the arclength derivative of x from the unit normal


% A 3d array containing the blocks for arclength differentiation
%ds = zeros(panel_res, panel_res, size(s, 2));
%for i = 1:size(s, 2)
%    ds(:,:,i) = (1./s(:,i)).*d_leg;
%end
% Construct a sparse block matrix for arclength differentiation
%i_func = @(i) panel_res.*floor(i./panel_res.^2) + mod(i, panel_res) + 1;
%j_func = @(i) floor(i./panel_res) + 1;
%i_vec = i_func(0:(numel(ds) - 1));
%j_vec = j_func(0:(numel(ds) - 1));
%ds_mat = sparse(i_vec, j_vec, ds(:));

%ds_vec = zeros(2*N);
%ds_vec(1:2:end, 1:2:end) = ds_mat;
%ds_vec(2:2:end, 2:2:end) = ds_mat;
%xs = ds_vec * chnkr_in.r(:);
%xss = ds_vec * xs(:); % Needed for the tension
%i = floor((2:(2*N + 1))./2);
%xs_mat = sparse(i, 1:(2*N), xs(:));
%xss_mat = sparse(i, 1:(2*N), xss(:));

k_st = STau_kernel(1);
k_s  = kernel("stokes", "s", 1);
S_mat = chunkermat(chnkr_in, k_s);
ST_mat = chunkermat(chnkr_in, k_st);

s_test = S_mat * chnkr_in.r(:);
st_test = ST_mat * chnkr_in.r(:);

%contract_1 = xs_mat * ST_mat * xs_mat';
%contract_2 = xs_mat * ST_mat * xss_mat';


A = CoefIntMat(panel_res);

tens_mat = make_tension_matrix(chnkr_in, ST_mat, smp, int_i);

tens_data = (reshape(chnkr_in.r(1, 1:end, 1:end), chnkr_in.npt, 1));
tens_int = (reshape(chnkr_in.r(2, 1:end, 1:end), chnkr_in.npt, 1));

for i = 1:pan_num
    tens_data((i - 1)*panel_res + 2:i*panel_res) = smp.dmat1 * tens_data((i - 1)*panel_res + 1:i*panel_res);
    tens_data((i - 1)*panel_res + 1) = smp.n.u(1, 1:end) * tens_int((i - 1)*panel_res + 1:i*panel_res);
end

bend_mat = make_bending_matrix(chnkr_in, ST_mat, smp, int_i);

%Test that x_s is a unit vector
%x_test = reshape(x_s, 2, N);
%test_norm = zeros(N, 1);
%for i = 1:N
%    test_norm(i) = norm(x_test(:, i)) - 1;
%end
[sys_mat, sys_dmat] = stokes_system_matrix(chnkr_in, ST_mat, smp, int_i, .1);
t_eval = tension_eval_matrix(chnkr_in, S_mat, smp, int_i);

[rho, L_output] = test_tension(chnkr_in, 4, smp);
[mu, B_output] = test_bending(chnkr_in, 4, smp, .1);

a = L_output - tens_mat * rho;
b = B_output - sys_mat * mu;

function [rho, L_output] = test_tension(chnkr, freq, smp)

    N = chnkr.npt;
    res = chnkr.k;
    p_num = N/res;

    theta = atan2(chnkr.r(2, 1:end), chnkr.r(1, 1:end))';
    rho = freq*cos(freq*theta);
    sigma = sin(freq*theta);
    for i = 1:p_num
        rho((i - 1)*res + 2:i*res) = smp.dmat1 * rho((i - 1)*res + 1:i*res);
        rho((i - 1)*res + 1) = smp.n.u(1, 1:end) * sigma((i - 1)*res + 1:i*res);
    end
    L_output = (-freq/4)*sigma;
    for i = 1:p_num
        L_output((i - 1)*res + 1:i*res) = smp.lob * smp.n.u * L_output((i - 1) * res + 1:i*res);
        L_output((i - 1)*res + 1) = 0;
    end
end

function [mu, B_output] = test_bending(chnkr, freq, smp, dt)
    N = chnkr.npt;
    res = chnkr.k;
    p_num = N/res;

    theta = atan2(chnkr.r(2, 1:end), chnkr.r(1, 1:end));
    x = [cos(freq * theta); sin(freq*theta)];
    x = x(:);
    x_s = [-freq*sin(freq * theta); freq * cos(freq*theta)];
    x_s = x_s(:);
    x_ss = [-freq^2*cos(freq*theta); -freq^2 * sin(freq*theta)];
    x_ss = x_ss(:);
    mu = [freq^3 * sin(freq*theta); -freq^3 * cos(freq*theta)];
    mu = mu(:);
    for i = 1:p_num
        j = 2*(i - 1)*res + 1;
        k = 2*i*res;
        mu(j + 6:2:k) = smp.dmat3 * mu(j:2:k);
        mu(j + 7:2:k) = smp.dmat3 * mu(j + 1:2:k);
        mu(j) = smp.n.u(1, 1:end) * x(j:2:k);
        mu(j + 1) = smp.n.u(1, 1:end) * x(j + 1:2:k);
        mu(j + 2) = smp.n.u(1, 1:end) * x_s(j:2:k);
        mu(j + 3) = smp.n.u(1, 1:end) * x_s(j + 1:2:k);
        mu(j + 4) = smp.n.u(1, 1:end) * x_ss(j:2:k);
        mu(j + 5) = smp.n.u(1, 1:end) * x_ss(j + 1:2:k);
    end
    B_output = (1 + .25*dt*freq^3)*x;
    for i = 1:p_num
        j = 2*(i - 1)*res + 1;
        k = 2*i*res;
        B_output(j + 4:2:k) = smp.l2.dmat * B_output(j:2:k);
        B_output(j + 5:2:k) = smp.l2.dmat * B_output(j + 1:2:k);
        B_output(j:j+5) = 0;
    end
end
