clear;
close all;
rng(42);

% ============== System parameters ============== 
Nt = 8;
nIRSrow = 8;
nIRScol = 8;
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
alphaMaxNum = 10.^((1:10)/10); % 1-10dB
% NsNumber = nIRSrow .* nIRScolNumber;

PropRISNumSimuRes = zeros(11, length(alphaMaxNum));
SDRRISNumSimuRes = zeros(11, length(alphaMaxNum));
MMRISNumSimuRes = zeros(11, length(alphaMaxNum));
PassiveRISNumSimuRes = zeros(11, length(alphaMaxNum));
PropStaticRISNumSimuRes = zeros(11, length(alphaMaxNum));

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

for n = 1 : length(alphaMaxNum)

    % nIRScol = nIRScolNumber(n);
    % nIRSrow = nIRScol;
    Ns = nIRSrow * nIRScol;
    alphaMax = alphaMaxNum(n);
    thetaVecCurrent = alphaMax * thetaVecCurrentDataSet{7} / 10;

    % thetaVecCurrent = thetaVecCurrentDataSet{7};

    fprintf(' ------- Current Number Of Alpha: %d ------- \n', alphaMax);

    resPropFinal = cell(100, 1);
    resSDRFinal = cell(100, 1);
    resMMFinal = cell(100, 1);
    resPropPassiveFinal = cell(100, 1);
    resPropStaticFinal = cell(100, 1);
    % --- Simulation ---
    tic
    for iter = 1 : 100
        fprintf(' === Monte Carlo Simulation: %d === \n', iter);
        % tic
        resPropFinal{iter} = mainFucProp(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, iter);
        % resPropStaticFinal{iter} = mainFucPropStatic(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, iter);
        % resSDRFinal{iter} = mainFucSDR(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, iter);
        % resMMFinal{iter} = mainFucMM(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, alphaMax, alpha, epsilon, thetaVecCurrent, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, iter);
        % resPropPassiveFinal{iter} = mainFucPropPassive(Nt, nIRSrow, nIRScol, Ns, K, R, Pt, P2, P, Lambda, sigma, 1, alpha, epsilon, thetaVecCurrent./alphaMax, G_NLoS{iter}, g_NLoS{iter}, h_NLoS{iter}, iter);
        % toc
    end
    toc
    PropStaticRISNumSimuRes(:, n) = mean(reshape(cell2mat(PropStaticFinal), 11, []), 2);
    PropRISNumSimuRes(:, n) = mean(reshape(cell2mat(resPropFinal), 11, []), 2);
    SDRRISNumSimuRes(:, n) = mean(reshape(cell2mat(resSDRFinal), 11, []), 2);
    MMRISNumSimuRes(:, n) = mean(reshape(cell2mat(resMMFinal), 11, []), 2);
    % PassiveRISNumSimuRes(:, n) = mean(reshape(cell2mat(resPropPassiveFinal), 21, []), 2);
end
% save("Prop_StaticRISNumRes.mat", "PropStaticRISNumSimuRes");
save("AlphaNumberRes.mat", "PropRISNumSimuRes", "SDRRISNumSimuRes", "MMRISNumSimuRes", "PropStaticRISNumSimuRes");