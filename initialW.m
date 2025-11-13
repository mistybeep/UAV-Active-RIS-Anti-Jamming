function wCurrent = initialW(Nt, K, dCurrent, gCurrent, h, G, g, thetaMatCurrent, Pt, P2, P)

    wCurrent = dCurrent * (dCurrent' * dCurrent)^(-1);
    % 
    wCurrent = sqrt(Pt) * wCurrent / norm(wCurrent);

end