function [peri, vol] = conserved_qs(chnkr)

    peri = sum(chnkr.wts(:));
    vol = sum(chnkr.wts(:) .* (chnkr.n(1, :) .* chnkr.r(1, :))');

end