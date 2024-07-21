function chnkr = update_chunker(chnkr, r_update)
    %method to update chunker geometry in case this is necessary
    % after modifying the boundary.
    % currently this does not maintain the parameterization, because
    % That requires calculating a quantity which is currently not
    % tracked by chunkers. 
    [~, ~, u, v] = lege.exps(chnkr.k);
    dermat2 = v*([lege.derpol(u); zeros(1, chnkr.k)]);
    chnkr.r(:) = r_update(:);
    for i = 1:(chnkr.npt/chnkr.k)
        chnkr.d(1:2, 1:end, i) = chnkr.r(1:2, 1:end, i) * dermat2';
        chnkr.d2(1:2, 1:end, i) = chnkr.d(1:2, 1:end, i) * dermat2';
    end
    chnkr.n(1, 1:end) = chnkr.d(2, 1:end);
    chnkr.n(2, 1:end) = -chnkr.d(1, 1:end);
    chnkr.n(1:2, :) = chnkr.n(1:2, :)./sqrt(sum(chnkr.d(1:2, :).^2, 1));
    chnkr.wts = weights(chnkr);
end