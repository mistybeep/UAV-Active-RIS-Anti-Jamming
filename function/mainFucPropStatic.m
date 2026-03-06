function tIter = mainFucPropStatic(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS, g_NLoS, h_NLoS, JamPosErr)
    
    % ============= initialization =============   
    iIter = 0;
    locCurrent = [50, 50, 60];

    % ============= IRS beamforming vector initialization  
    thetaMatCurrent = diag(thetaVecCurrent);
    
    % ============= Channel generation and normalization 
    [G, GBar, g, gBar, h, hBar, locU, locJam, locBS] = ChanGen(Nt, K, R, nIRSrow, nIRScol, locCurrent, G_NLoS, g_NLoS, h_NLoS, Lambda);
    G = G / sigma;
    g = g / sigma;
    dCurrent = (h' * thetaMatCurrent * G)';
    gCurrent = h' * thetaMatCurrent * g;
    % ============= Transmit beamformer initialization
    wCurrent = initialW(Nt, K, dCurrent, gCurrent, h, G, g, thetaMatCurrent, Pt, P2, P);
    [tCurrent, muCurrent] = K_update_SINR(wCurrent, dCurrent, gCurrent, thetaMatCurrent, h, P, K);
    locJamTrue = locJam + JamPosErr;
    % ============= 
    realChange = 1;
    
    tIter = zeros(11, 1);
    tTemp = zeros(11, 1);
    
    tIter(1) = 10*log10(tCurrent);
    tTemp(1) = 10*log10(tCurrent);
    
    while realChange >= epsilon && iIter < 10
    
        iIter = iIter+1;
        % 
        [G, GBar, g, gBar, h, hBar, locU, locJam, locBS] = ChanGen(Nt, K, R, nIRSrow, nIRScol, locCurrent, G_NLoS, g_NLoS, h_NLoS, Lambda);
        G = G / sigma;
        g = g / sigma;
        dCurrent = (h' * thetaMatCurrent * G)';
        gCurrent = h' * thetaMatCurrent * g;

        IterD = 0;
        while IterD < 20
            IterD = IterD + 1;
            [W, thetaMat] = updateBeam(Nt, K, Ns, R, thetaMatCurrent, wCurrent, dCurrent, tCurrent, muCurrent, G, g, h, Pt, P2, P, alphaMax);
    
            wCurrent = W;
            thetaMatCurrent = thetaMat;
            dCurrent = (h' * thetaMatCurrent * G)';
            gCurrent = h' * thetaMatCurrent * g;
            [tNew, muNew] = K_update_SINR(wCurrent, dCurrent, gCurrent, thetaMatCurrent, h, P, K);
            if abs(10*log10(tNew) - 10*log10(tCurrent)) > 1e-3
                tCurrent = tNew;
                muCurrent = muNew;
            else
                break;
            end
        end

        realChange = abs(10*log10(tCurrent) - tIter(iIter));
    
        [G, GBar, g, gBar, h, hBar, locU, ~, locBS] = ChanGen2(locJamTrue, Nt, K, R, nIRSrow, nIRScol, locCurrent, G_NLoS, g_NLoS, h_NLoS, Lambda);
        G = G / sigma;
        g = g / sigma;
        dCurrent = (h' * thetaMatCurrent * G)';
        gCurrent = h' * thetaMatCurrent * g;
        [SINRTemp, ~] = K_update_SINR(wCurrent, dCurrent, gCurrent, thetaMatCurrent, h, P, K);
        tTemp(iIter + 1) = 10*log10(SINRTemp);
        tIter(iIter + 1) = 10*log10(tCurrent);

    end
    tTemp(iIter + 1:end) = 10*log10(SINRTemp);
    tIter = tTemp;

end