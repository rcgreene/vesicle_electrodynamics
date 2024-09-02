function chnkr = update_chunker(varargin)
    %method to update chunker geometry in case this is necessary
    % after modifying the boundary.
    % currently this does not maintain the parameterization, because
    % That requires calculating a quantity which is currently not
    % tracked by chunkers. 

    % first argument is the chnkr to modify, second is the new chnkr.r.
    % If no second argument is passed, then d, d2, n, and wts will be
    % redetermined anyway. This is for instances where chnkr.r is changed
    % outside of the routine.

    chnkr = varargin{1};
    if nargin == 2
        chnkr.r = varargin{2};
    end
    [~, ~, u, v] = lege.exps(chnkr.k);
    dermat2 = v*([lege.derpol(u); zeros(1, chnkr.k)]);
    for i = 1:(chnkr.npt/chnkr.k)
        chnkr.d(1:2, 1:end, i) = chnkr.r(1:2, 1:end, i) * dermat2';
        chnkr.d2(1:2, 1:end, i) = chnkr.d(1:2, 1:end, i) * dermat2';
    end
    chnkr.n(1, 1:end) = chnkr.d(2, 1:end);
    chnkr.n(2, 1:end) = -chnkr.d(1, 1:end);
    chnkr.n(1:2, :) = chnkr.n(1:2, :)./sqrt(sum(chnkr.d(1:2, :).^2, 1));
    chnkr.wts = weights(chnkr);
end