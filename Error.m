clear;
close all;
rng(42);

% ============== System parameters ============== 
Nt = 8;
nIRSrow = 4;
nIRScol = 4;
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

% --- 并行设置 ---
% delete(gcp('nocreate'));
% numWorkers = 15;
% parpool('local', numWorkers);
addpath("function\");
error = 0:2:20;
PropRISNumSimuRes = zeros(11, length(error));
SDRRISNumSimuRes = zeros(11, length(error));
MMRISNumSimuRes = zeros(11, length(error));
PassiveRISNumSimuRes = zeros(11, length(error));
PropStaticRISNumSimuRes = zeros(11, length(error));

load('RISNumDataSet.mat', "thetaVecCurrentDataSet", "G_NLoSDataSet", "g_NLoSDataSet", "h_NLoSDataSet");

thetaVecCurrent = thetaVecCurrentDataSet{3};
G_NLoS = cell(100, 1);
g_NLoS = cell(100, 1);
h_NLoS = cell(100, 1);
JamErr = cell(100, 1);
JamErrRad = cell(100, 1);
for i = 1 : 100
    G_NLoS{i} = G_NLoSDataSet{i, 3};
    g_NLoS{i} = g_NLoSDataSet{i, 3};
    h_NLoS{i} = h_NLoSDataSet{i, 3};
    error_dir = randn(R, 2);
    JamErr{i} = [error_dir, zeros(R, 1)];
    JamErrRad{i} = rand();
end


for n = 9 : length(error)

    % nIRScol = nIRScolNumber(n);
    % nIRSrow = nIRScol;
    Ns = nIRSrow * nIRScol;

    thetaVecCurrent = thetaVecCurrentDataSet{3};
    Err = error(n);
    fprintf(' ------- Current Number Of Error: %d ------- \n', Err);

    resPropFinal = cell(100, 1);
    resSDRFinal = cell(100, 1);
    resMMFinal = cell(100, 1);
    resPropPassiveFinal = cell(100, 1);
    resPropStaticFinal = cell(100, 1);
    
    % locJam_est = locJam + error_vec; % 系统已知的“估计干扰位置”
    % JamPosErr = zeros(R, 3);
    
    % --- Simulation ---
    tic
    parfor iter = 1 : 100
        
        fprintf(' === Monte Carlo Simulation: %d === \n', iter);
        err = (JamErr{iter} / norm(JamErr{iter})) * Err * JamErrRad{iter};
        % tic
        resPropFinal{iter} = mainFucProp(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, err);
        resPropStaticFinal{iter} = mainFucPropStatic(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, err);
        resSDRFinal{iter} = mainFucSDR(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, err);
        resMMFinal{iter} = mainFucMM(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, err);
        resPropPassiveFinal{iter} = mainFucPropPassive(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, 1, alpha, epsilon, thetaVecCurrent./alphaMax, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, err);
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
save("ErrorNumberRes.mat", "PropRISNumSimuRes", "SDRRISNumSimuRes", "MMRISNumSimuRes", "PropStaticRISNumSimuRes", "PassiveRISNumSimuRes");