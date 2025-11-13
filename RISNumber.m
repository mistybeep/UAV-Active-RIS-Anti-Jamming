clear;
close all;
rng(42);

% ============== System parameters ============== 
Nt = 8;
nIRSrow = 4;
nIRScol = 9;
Ns = nIRSrow*nIRScol;
K = 3;
R = 2;
Pt = db2pow(40 - 30);
P2 = db2pow(50 - 30);
P = 500; % -20dB
f = 3e9;
c = 3e8;
Lambda = c/f;
N0 = db2pow(-174-30);
B = 20e6;
% sigma = sqrt(B*N0);
sigma = 1e-5; % -70dBm
alphaMax = 10; % amplification gain
alpha = 2.2;
epsilon = 1e-2;
relChange = 1e3;

% --- 并行设置 ---
delete(gcp('nocreate'));
numWorkers = 15;
parpool('local', numWorkers);
addAttachedFiles(gcp, {'D:\Software\MATLAB\R2023b\bin\cvx'});

nIRScolNumber = 4:9;
% NsNumber = nIRSrow .* nIRScolNumber;

PropRISNumSimuRes = zeros(length(nIRScolNumber), 1);

load('RISNumDataSet.mat', "thetaVecCurrentDataSet", "G_NLoSDataSet", "g_NLoSDataSet", "h_NLoSDataSet");

for n = 1 : length(nIRScolNumber)

    nIRScol = nIRScolNumber(n);
    Ns = nIRSrow * nIRScol;

    thetaVecCurrent = thetaVecCurrentDataSet{n};

    fprintf(' ------- Current Number Of Active RIS: %d ------- \n', Ns);

    resPropFinal = cell(100, 1);
    % --- Simulation ---
    tic
    parfor iter = 1 : 100
        fprintf(' === Monte Carlo Simulation: %d === \n', iter);
        resPropFinal{iter} = mainFucProp(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoSDataSet{iter, n}, g_NLoSDataSet{iter, n}, h_NLoSDataSet{iter, n}, iter);
    end
    toc
    PropRISNumSimuRes(n) = mean(cell2mat(resPropFinal));
end

save("RISNumberRes.mat", "PropRISNumSimuRes");

function res = mainFucProp(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS, g_NLoS, h_NLoS, iter)
    
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
    
    % ============= 
    realChange = 1;
    
    tIter = 10*log10(tCurrent);
    
    while realChange >= epsilon && iIter < 20
    
        iIter = iIter+1;
        
        IterD = 0;
        while IterD < 20
            IterD = IterD + 1;
            [locCurrent, G, g, h] = updateDeployment(Nt, K, Ns, R, GBar, gBar, hBar, wCurrent, thetaMatCurrent, tCurrent, muCurrent, locCurrent', locU', locJam', locBS', alpha, P, P2);
            locCurrent = locCurrent';
            % [G, GBar, g, gBar, h, hBar, locU, locJam, locBS] = ChanGen(Nt, K, R, nIRSrow, nIRScol, locCurrent, G_NLoS, g_NLoS, h_NLoS, Lambda);
            G = G / sigma;
            g = g / sigma;
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

        realChange = abs(10*log10(tCurrent) - tIter(end));
    
        tIter = [tIter, 10*log10(tCurrent)];

    end

    [tCurrent, ~] = K_update_SINR(wCurrent, dCurrent, gCurrent, thetaMatCurrent, h, P, K);
    res = 10*log10(tCurrent);

end