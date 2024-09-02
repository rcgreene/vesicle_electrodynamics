function [output_mat, shift_mat] = tens_solve(chnkr, stt_kernel_mat, smp, int_i)

panel_res = chnkr.k;
pan_num = chnkr.npt/panel_res;
N = chnkr.npt;


i = 1:2*N;
j = floor((i + 1)./2);
xsl_mat = sparse(j, i, int_i.l.xs);

% resample a panel at lobatto nodes
tmp = zeros(2*panel_res, 2*panel_res);
tmp(1:2:end, 1:2:end) = smp.lob *  smp.n.u;
tmp(2:2:end, 2:2:end) = smp.lob *  smp.n.u;
% resample for each panel
lob_resample = kron(speye(pan_num), tmp);


% inverse resample. Lobatto to legendre
%lob_i_resample = kron(speye(pan_num), inv(tmp));

% Downsample to n - 1 legendre nodes
%tmp(3:2:end, 1:2:end) = smp.n_1.v * eye(panel_res - 1, panel_res)*smp.n.u;
%tmp(4:2:end, 2:2:end) = smp.n_1.v * eye(panel_res - 1, panel_res)*smp.n.u;
%tmp(1:2, 1:end) = 0;
%dsamp_mat = kron(speye(pan_num), tmp);

% Block diagonal N x 2N matrices lifting a projecting a vector field to a
% scalar
xs_mat = sparse(j, i, int_i.n.xs);
xss_mat = sparse(j, i, int_i.n.xss);

%xs_n_1 = zeros(2, panel_res, pan_num);
%xs_n_1(1:2, 2:end, 1:end) = reshape(int_i.n_1.xs, 2, panel_res - 1, pan_num);
%xs_n_1_mat = sparse(j, i, xs_n_1(:));
shift_mat = 0; %temporary   %xs_n_1_mat * dsamp_mat * stt_kernel_mat' * xs_mat';
c_1 = xsl_mat * lob_resample * stt_kernel_mat * xs_mat';
c_2 = xsl_mat * lob_resample * stt_kernel_mat * xss_mat';

%c_2 = resample * shift_mat * xs_mat * stt_kernel_mat * xss_mat';
%c_1 = resample * shift_mat * shift_mat;
%shift_mat = shift_mat;

us_mat = zeros(panel_res);
us_mat(1:end, 2:end) = smp.umat1;
output_mat = c_1 * kron(speye(pan_num), us_mat);
a = 0:((panel_res^2)*pan_num - 1);
i = (panel_res)*floor(a/(panel_res^2)) + mod(a, panel_res) + 1;
j = (panel_res)*floor(a/(panel_res^2)) + floor(mod(a, panel_res.^2)/panel_res) + 1;

output_mat = output_mat + ...
    c_2 * kron(speye(pan_num), smp.n.v) * sparse(i, j, int_i.int1(:));
%tmp = zeros(panel_res, panel_res);
%tmp(2:end, 1:end) = smp.dmat1;
%output_mat = kron(speye(pan_num), tmp) * output_mat;
%for i = 1:pan_num
%    j = chnkr.adj(2, i) - 1;
%    output_mat(i*panel_res, 1:end) = ...
%        (output_mat(i*panel_res, 1:end) + output_mat(j*panel_res + 1, 1:end))/2;
%end
for i = 0:pan_num-1
    output_mat(i*panel_res + 1, 1:end) = 0;
    if i == 0
        output_mat(1, 1:end) = reshape(int_i.int1(2, 1:end, 1:end), 1, N);
    else
        %j = chnkr.adj(2, i + 1);
        %output_mat(i*panel_res + 1, 1:end) = ...
        %    output_mat(i * panel_res + 1, 1:end) ...
        %    -output_mat(j*panel_res, 1:end);
        output_mat(i*panel_res + 1, i*panel_res + 1:(i + 1)*panel_res) = ...
            sum(int_i.int1(1:end, 1:end, i + 1), 1);
        j = chnkr.adj(2, i + 1) - 1;
        output_mat(i*panel_res + 1, j*panel_res + 1:(j + 1)*panel_res) = ...
            -sum(int_i.int1(1:2:end, 1:end, j + 1), 1) + sum(int_i.int1(2:2:end, 1:end, j + 1), 1);
    end
end

end