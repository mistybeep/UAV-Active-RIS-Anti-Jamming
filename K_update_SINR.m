function [SINR, num] = K_update_SINR(W, dCurrent, gCurrent, thetaMatCurrent, h, P, K)

    He = abs(dCurrent' * W).^2;
    gr = zeros(K, 1);
    num = zeros(K, 1);
    for k = 1 : K
        tmp=He(k,:);
        num(k) = (sum(tmp)-tmp(k) + P * sum(abs(gCurrent(k, :)).^2) + norm(h(:, k)' * thetaMatCurrent, 2)^2 + 1);
        gr(k)=tmp(k)/num(k);
    end

    SINR = min(gr);
end