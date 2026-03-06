clear;
close all;
rng(42);

% ============== System parameters
Nt = 8;
nIRSrow = 4;
nIRScol = 4;
Ns = nIRSrow*nIRScol;
K = 3;
R = 2;
Pt = db2pow(30 - 30); % 30dBm
P2 = db2pow(20 - 30); % 20dBm
f = 3e9;
c = 3e8;
Lambda = c/f;
sigma = 1e-5;
alphaMax = 10; % amplification gain
alpha = 2.2;
epsilon = 1e-2;
relChange = 1e3;
addpath("function\");
% --- 并行设置 ---
% delete(gcp('nocreate'));
% numWorkers = 12;
% parpool('local', numWorkers);
% addAttachedFiles(gcp, {'D:\Software\MATLAB\R2023b\bin\cvx'});

PNum = Pt./(R * [1e-2, 1e-1, 1]);
PropJamPowSimuRes = zeros(11, length(PNum));
SDRJamPowSimuRes = zeros(11, length(PNum));
MMJamPowSimuRes = zeros(11, length(PNum));
PassiveJamPowSimuRes = zeros(11, length(PNum));
PropStaticJamPowSimuRes = zeros(11, length(PNum));

load('RISNumDataSet.mat', "thetaVecCurrentDataSet", "G_NLoSDataSet", "g_NLoSDataSet", "h_NLoSDataSet");

thetaVecCurrent = thetaVecCurrentDataSet{7};
G_NLoS = cell(100, 1);
g_NLoS = cell(100, 1);
h_NLoS = cell(100, 1);
for i = 1 : 100
    G_NLoS{i} = G_NLoSDataSet{i, 7};
    g_NLoS{i} = g_NLoSDataSet{i, 7};
    h_NLoS{i} = h_NLoSDataSet{i, 7};
end

for n = 1 : length(PNum)

    P = PNum(n);

    fprintf(' ------- Current Jammer Power: %d W ------- \n', P);

    resPropFinal = cell(100, 1);
    resSDRFinal = cell(100, 1);
    resMMFinal = cell(100, 1);
    resPropPassiveFinal = cell(100, 1);
    resPropStaticFinal = cell(100, 1);
    % --- Simulation ---
    tic
    parfor iter = 1 : 100
        fprintf(' === Monte Carlo Simulation: %d === \n', iter);
        % tic
        resPropFinal{iter} = mainFucProp(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, 0);
        resPropStaticFinal{iter} = mainFucPropStatic(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, 0);
        resSDRFinal{iter} = mainFucSDR(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, 0);
        resMMFinal{iter} = mainFucMM(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, 0);
        resPropPassiveFinal{iter} = mainFucPropPassive(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, 1, alpha, epsilon, thetaVecCurrent./alphaMax, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, 0);
        % toc
    end
    toc
    PropStaticJamPowSimuRes(:, n) = mean(reshape(cell2mat(resPropStaticFinal), 11, []), 2);
    PropJamPowSimuRes(:, n) = mean(reshape(cell2mat(resPropFinal), 11, []), 2);
    SDRJamPowSimuRes(:, n) = mean(reshape(cell2mat(resSDRFinal), 11, []), 2);
    MMJamPowSimuRes(:, n) = mean(reshape(cell2mat(resMMFinal), 11, []), 2);
    PassiveJamPowSimuRes(:, n) = mean(reshape(cell2mat(resPropPassiveFinal), 11, []), 2);
end

save("JamPowerRes2.mat", "PropJamPowSimuRes", "SDRJamPowSimuRes", "MMJamPowSimuRes", "PassiveJamPowSimuRes", "PropStaticJamPowSimuRes");
% save("Prop_StaticJamPowNumRes.mat", "PropStaticJamPowerNumSimuRes");