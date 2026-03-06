clear;
close all;
rng(2026);

% ============== System parameters
Nt = 8;
nIRSrow = 4;
nIRScol = 9;
Ns = nIRSrow*nIRScol;
K = 3;
R = 2;
Pt = db2pow(30 - 30); % 30dBm
P2 = db2pow(20 - 30); % 40dBm
P = 500; % 40dBm
f = 3e9;
c = 3e8;
Lambda = c/f;
N0 = db2pow(-174-30);
B = 20e6;
% sigma = sqrt(B*N0);
sigma = 1e-5; % -70dBm
alphaMax = 10; % amplification gain 10dB
alpha = 2.2;
epsilon = 1e-2;
relChange = 1e3;

% --- 并行设置 ---
% delete(gcp('nocreate'));
% numWorkers = 10;
% parpool('local', numWorkers);
% addAttachedFiles(gcp, {'D:\Software\MATLAB\R2023b\bin\cvx'});

PNum = Pt./(R * [1, 1e-2]); % 40dBm
PropStaticJamPowerNumSimuRes = zeros(length(PNum), 1);

load('RISNumDataSet.mat', "thetaVecCurrentDataSet", "G_NLoSDataSet", "g_NLoSDataSet", "h_NLoSDataSet");

thetaVecCurrent = thetaVecCurrentDataSet{6};
G_NLoS = cell(100, 1);
g_NLoS = cell(100, 1);
h_NLoS = cell(100, 1);
for i = 1 : 100
    G_NLoS{i} = G_NLoSDataSet{i, 6};
    g_NLoS{i} = g_NLoSDataSet{i, 6};
    h_NLoS{i} = h_NLoSDataSet{i, 6};
end

nIRScolNumber = 2:9;
% ============= NLOS ===============
G_NLoSDataSet = cell(100, length(nIRScolNumber));
g_NLoSDataSet = cell(100, length(nIRScolNumber));
h_NLoSDataSet = cell(100, length(nIRScolNumber));

% ============= RIS Init ===============
thetaVecCurrentDataSet = cell(length(nIRScolNumber), 1);

for n = 1 : length(nIRScolNumber)

    % nIRSrow = nIRSrowNumber(n);
    nIRScol = nIRScolNumber(n);
    nIRSrow = nIRScol;
    Ns = nIRSrow * nIRScol;

    thetaVecCurrentDataSet{n} = alphaMax * exp(1j*2*pi*rand(Ns, 1));

    for i = 1 : 100
        G_NLoSDataSet{i, n} = sqrt(1/2) * (randn(nIRSrow*nIRScol, Nt) + 1j*randn(nIRSrow*nIRScol, Nt));
        g_NLoSDataSet{i, n} = sqrt(1/2) * (randn(nIRSrow*nIRScol, R) + 1j*randn(nIRSrow*nIRScol, R));
        h_NLoSDataSet{i, n} = sqrt(1/2) * (randn(K, nIRSrow*nIRScol) + 1j*randn(K, nIRSrow*nIRScol))';
    end
end

save('RISNumDataSet.mat', "thetaVecCurrentDataSet", "G_NLoSDataSet", "g_NLoSDataSet", "h_NLoSDataSet");

% load('RISNumDataSet.mat', "thetaVecCurrentDataSet", "G_NLoSDataSet", "g_NLoSDataSet", "h_NLoSDataSet");

for n = 1 : length(PNum)

    P = PNum(n);

    fprintf(' ------- Current Jammer Power: %d W ------- \n', P);
    resFinal = cell(100, 1);
    % --- Simulation ---
    tic
    for iter = 1 : 100
        
        resFinal{iter} = mainFuc(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, iter);
        
    end
    toc
    PropStaticJamPowerNumSimuRes(n) = mean(cell2mat(resFinal));
end

save("Prop_StaticJamPowNumRes.mat", "PropStaticJamPowerNumSimuRes");

function [res, tIter] = mainFuc(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS, g_NLoS, h_NLoS, iter)
    
    % ============= initialization =============   
    iIter = 0;
    locCurrent = [50, 50, 60];

    fprintf(' === Monte Carlo Simulation: %d === \n', iter);

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
        while IterD < 10
            IterD = IterD + 1;
            [locCurrent, G, g, h] = updateDeployment(Nt, K, Ns, R, GBar, gBar, hBar, wCurrent, thetaMatCurrent, tCurrent, muCurrent, locCurrent', locU', locJam', locBS', alpha, P, P2);
            locCurrent = locCurrent';
            % [G, GBar, g, gBar, h, hBar, locU, locJam, locBS] = ChanGen(Nt, K, R, nIRSrow, nIRScol, locCurrent, G_NLoS, g_NLoS, h_NLoS, Lambda);
            G = G / sigma;
            g = g / sigma;
            dCurrent = (h' * thetaMatCurrent * G)';
            gCurrent = h' * thetaMatCurrent * g;

            [tNew, muNew] = K_update_SINR(wCurrent, dCurrent, gCurrent, thetaMatCurrent, h, P, K);
            if abs(tNew - tCurrent) > 1e-3
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
        while IterD < 10
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

        [tCurrent, ~] = K_update_SINR(wCurrent, dCurrent, gCurrent, thetaMatCurrent, h, P, K);

        realChange = abs(10*log10(tCurrent) - tIter(end));
    
        tIter = [tIter, 10*log10(tCurrent)];

    end
    
    [tCurrent, ~] = K_update_SINR(wCurrent, dCurrent, gCurrent, thetaMatCurrent, h, P, K);
    res = 10*log10(tCurrent);

end