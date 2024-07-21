function output_mat = tens_mat(chnkr, stt_kernel_mat, smp, int_i)

panel_res = chnkr.k;
pan_num = chnkr.npt/panel_res;
N = chnkr.npt;

tens_mat = zeros(N, N);

i = 1:2*N;
j = floor((i + 1)./2);
xsl_mat = sparse(j, i, int_i.l.xs);
lob_resample = kron(speye(pan_num), smp.lob * smp.n.u);
xs_mat = sparse(j, i, int_i.n.xs);
xss_mat = sparse(j, i, int_i.n.xss);
c_1 = xsl_mat * lob_resample * stt_kernel_mat * xs_mat';
c_2 = xsl_mat * lob_resample * stt_kernel_mat * xss_mat';

us_mat = zeros(panel_res);
us_mat(1:end, 2:end) = smp.umat1;
end