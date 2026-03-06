clear;
close all;
rng(2026);

% ============== System parameters ============== 
Nt = 8;
nIRSrow = 4;
nIRScol = 9;
Ns = nIRSrow*nIRScol;
K = 3;
R = 2;
Pt = db2pow(30 - 30);
P2 = db2pow(20 - 30);
P = Pt./(R * 1e-2);
f = 3e9;
c = 3e8;
Lambda = c/f;
sigma = 1e-5; % -70dBm
alphaMax = 10; % amplification gain
alpha = 2.2;
epsilon = 1e-2;
relChange = 1e3;

addpath("function\");
% --- 并行设置 ---
% delete(gcp('nocreate'));
% numWorkers = 12;
% parpool('local', numWorkers);

nIRScolNumber = 2:9;
% NsNumber = nIRSrow .* nIRScolNumber;

PropRISNumSimuRes = zeros(11, length(nIRScolNumber));
SDRRISNumSimuRes = zeros(11, length(nIRScolNumber));
MMRISNumSimuRes = zeros(11, length(nIRScolNumber));
PassiveRISNumSimuRes = zeros(11, length(nIRScolNumber));
PropStaticRISNumSimuRes = zeros(11, length(nIRScolNumber));

load('RISNumDataSet.mat', "thetaVecCurrentDataSet", "G_NLoSDataSet", "g_NLoSDataSet", "h_NLoSDataSet");

for n = 1 : length(nIRScolNumber)

    nIRScol = nIRScolNumber(n);
    nIRSrow = nIRScol;
    Ns = nIRSrow * nIRScol;

    thetaVecCurrent = thetaVecCurrentDataSet{n};

    fprintf(' ------- Current Number Of Active RIS: %d ------- \n', Ns);

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
        resPropFinal{iter} = mainFucProp(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoSDataSet{iter, n}, g_NLoSDataSet{iter, n}, h_NLoSDataSet{iter, n}, iter);
        % % 
        resPropStaticFinal{iter} = mainFucPropStatic(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoSDataSet{iter, n}, g_NLoSDataSet{iter, n}, h_NLoSDataSet{iter, n}, iter);

        resSDRFinal{iter} = mainFucSDR(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoSDataSet{iter, n}, g_NLoSDataSet{iter, n}, h_NLoSDataSet{iter, n}, iter);
        % 
        resMMFinal{iter} = mainFucMM(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoSDataSet{iter, n}, g_NLoSDataSet{iter, n}, h_NLoSDataSet{iter, n}, iter);

        resPropPassiveFinal{iter} = mainFucPropPassive(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, 1, alpha, epsilon, thetaVecCurrent./alphaMax, G_NLoSDataSet{iter, n}, g_NLoSDataSet{iter, n}, h_NLoSDataSet{iter, n}, iter);
        % toc
    end
    toc
    PropStaticRISNumSimuRes(:, n) = mean(reshape(cell2mat(resPropStaticFinal), 11, []), 2);
    PropRISNumSimuRes(:, n) = mean(reshape(cell2mat(resPropFinal), 11, []), 2);
    SDRRISNumSimuRes(:, n) = mean(reshape(cell2mat(resSDRFinal), 11, []), 2);
    MMRISNumSimuRes(:, n) = mean(reshape(cell2mat(resMMFinal), 11, []), 2);
    PassiveRISNumSimuRes(:, n) = mean(reshape(cell2mat(resPropPassiveFinal), 11, []), 2);
end
% save("Prop_StaticRISNumRes.mat", "PropStaticRISNumSimuRes");
save("RISNumberRes.mat", "PropRISNumSimuRes", "SDRRISNumSimuRes", "MMRISNumSimuRes", "PassiveRISNumSimuRes", "PropStaticRISNumSimuRes");