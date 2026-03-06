function [G, GBar, g, gBar, h, hBar, locU, locJam, locBS] = ChanGen2(locJam, Nt, K, R, nIRSrow, nIRScol, locUAV, G_NLoS, g_NLoS, h_NLoS, Lambda)
    halfLambda = 0.5*Lambda;
    quarterLambda = 0.25*Lambda;
    kappa = db2pow(3);
    alpha = 2.2;
    beta = db2pow(-30.18);
    %=========== Location of nodes/antennas/tiles (all in m)
    %----------- tx uniform linear array (ULA) G_NLoS, g_NLoS, h_NLoS, 
    locBS = zeros(1, 3);
    xt = 50;
    yt = 0;
    zt = 10;
    locBS(1) = xt;
    locBS(2) = yt;
    locBS(3) = zt;

    %----------- user single antenna
    locU = zeros(K, 3);
    locU(1, 1) = 25;
    locU(1, 2) = 60;
    locU(1, 3) = 2; 

    locU(2, 1) = 50;
    locU(2, 2) = 75;
    locU(2, 3) = 2; 

    locU(3, 1) = 75;
    locU(3, 2) = 60;
    locU(3, 3) = 2; 

    %----------- IRS uniform planar array (UPA)
    xs = locUAV(1);
    ys = locUAV(2);
    zs = locUAV(3);

    %----------- jammer single antenna
    % locJam = zeros(R, 3);
    % locJam(1, 1) = 30;
    % locJam(1, 2) = 30;
    % locJam(1, 3) = 2; 
    % 
    % locJam(2, 1) = 70;
    % locJam(2, 2) = 30;
    % locJam(2, 3) = 2; 

    %--- Transmitter (Tx) antenna array coordinates ---
    % locTcenter = [xt, yt, zt];  % 1x3 vector
    % locT = repmat(locTcenter, Nt, 1);  % Nt x 3
    % 
    % if mod(Nt, 2) == 0
    %     locT(1, 3) = 0 - 0.5 * (Nt - 2) * halfLambda - quarterLambda;
    % else
    %     locT(1, 3) = 0 - 0.5 * (Nt - 1) * halfLambda;
    % end
    % 
    % % Generate y-coordinates for all Nt antennas
    % z_vals = locT(1, 3) + (0:Nt-1)' * halfLambda;  % column vector
    % locT(:, 3) = z_vals;
    locT = zeros(Nt, 1);
    if mod(Nt, 2) == 0
        locT(1) = 0 - 0.5 * (Nt - 2) * halfLambda - quarterLambda;
    else
        locT(1) = 0 - 0.5 * (Nt - 1) * halfLambda;
    end
    z_vals = locT(1) + (0:Nt-1)' * halfLambda;
    locT = z_vals;

    %--- IRS (Intelligent Reflecting Surface) coordinates ---
    % locIRScenter = locUAV;  % 1x3
    % Initialize 3D array: nIRSrow x nIRScol x 3
    % locS = zeros(nIRSrow, nIRScol, 3);
    % locS(:, :, 1) = locUAV(1);
    % locS(:, :, 2) = locUAV(2);
    % locS(:, :, 3) = locUAV(3);
    % 
    % if mod(nIRScol, 2) == 0
    %     locS(:, :, 1) = xs - 0.5 * (nIRScol - 2) * halfLambda - quarterLambda;
    % else
    %     locS(:, :, 1) = xs - 0.5 * (nIRScol - 1) * halfLambda;
    % end
    % 
    % % Now assign x-coordinates for each column
    % x_start = locS(1, 1, 1);  % starting x for first column
    % x_vals = x_start + (0:nIRScol-1) * halfLambda;  % 1 x nIRScol
    % 
    % for nRow = 1:nIRSrow
    %     locS(nRow, :, 1) = x_vals;
    % end
    % 
    % % --- Set z-coordinates (rows, along dimension 3) ---
    % if mod(nIRSrow, 2) == 0
    %     y_start = ys - 0.5 * (nIRSrow - 2) * halfLambda - quarterLambda;
    % else
    %     y_start = ys - 0.5 * (nIRSrow - 1) * halfLambda;
    % end
    % y_vals = y_start + (0:nIRSrow-1)' * halfLambda;  % nIRSrow x 1
    % 
    % for jCol = 1:nIRScol
    %     locS(:, jCol, 2) = y_vals;
    % end
    locS = zeros(nIRSrow, nIRScol, 2);
    if mod(nIRScol, 2) == 0
        x_start = 0 - 0.5 * (nIRScol - 2) * halfLambda - quarterLambda;
    else
        x_start = 0 - 0.5 * (nIRScol - 1) * halfLambda;
    end

    x_vals = x_start + (0:nIRScol-1) * halfLambda;  % 1 x nIRScol

    for nRow = 1:nIRSrow
        locS(nRow, :, 1) = x_vals;
    end
    if mod(nIRSrow, 2) == 0
        y_start = 0 - 0.5 * (nIRSrow - 2) * halfLambda - quarterLambda;
    else
        y_start = 0 - 0.5 * (nIRSrow - 1) * halfLambda;
    end
    y_vals = y_start + (0:nIRSrow-1)' * halfLambda;  % nIRSrow x 1

    for jCol = 1:nIRScol
        locS(:, jCol, 2) = y_vals;
    end

    % --- Reshape IRS coordinates to (nIRSrow * nIRScol) x 3 matrix ---
    locS = reshape(locS, [nIRSrow * nIRScol, 2]);
    
    % --- dSU: distance from each user (K) to IRS -> K x 1
    dSU = sqrt(sum((reshape(locU, K, 1, 3) - reshape(locUAV, 1, 1, 3)).^2, 3));
    % % 
    % % --- dTS: distance from each Tx antenna (Nt) to each IRS element (Ns) -> Nt x Ns
    % dTS = sqrt(sum((reshape(locT, Nt, 1, 3) - reshape(locS, 1, nIRSrow * nIRScol, 3)).^2, 3))';
    % 
    % --- dJS: distance from Jammer to each IRS element (Ns) -> Nt x Ns
    dJS = sqrt(sum((reshape(locJam, R, 1, 3) - reshape(locUAV, 1, 1, 3)).^2, 3));

    % --- Tx to RIS channel
    theta = asin((locUAV(3)-zt)/norm(locUAV - [xt, yt, zt]));
    phi = atan2(locUAV(2) - yt, locUAV(1) - xt);
    
    G_Los = (exp(1j*2*pi*(cos(theta)*cos(phi) * locS(:, 1) + cos(theta)*sin(phi) * locS(:, 2) ) / Lambda))*(exp(1j*2*pi*(sin(theta) * locT) / Lambda))';
    % G_Los = exp(1j*2*pi*dTS/Lambda);
    GBar = sqrt(beta) * (sqrt(kappa/(kappa+1))*G_Los + sqrt(1/(kappa+1))*G_NLoS);
    G = sqrt(norm(locUAV - [xt, yt, zt])^(-alpha)) *  GBar;

    % --- Jam to RIS channel
    theta = asin((locUAV(3)*ones(R, 1)-locJam(:,3))./dJS);
    phi = atan2(locUAV(2)*ones(R, 1) - locJam(:, 2), locUAV(1)*ones(R, 1) - locJam(:, 1));

    % g_Los = exp(1j*2*pi*dJS/Lambda);
    g_Los = exp(1j*2*pi*(locS(:, 1) * (cos(theta).*cos(phi))' + locS(:, 2)* (cos(theta).*sin(phi))') / Lambda);
    % exp(1j*2*pi*(99.9338)/Lambda)
    gBar = sqrt(beta) * (sqrt(kappa/(kappa+1))*g_Los + sqrt(1/(kappa+1))*g_NLoS);
    g = gBar .* repmat(sqrt(sqrt(sum((repmat(locUAV, R, 1) - locJam).^2, 2)).^(-alpha))', nIRSrow * nIRScol, 1);

    % --- RIS to user channel
    theta = asin(-(locUAV(3)*ones(K, 1)-locU(:,3))./dSU);
    phi = atan2(locU(:, 2) - locUAV(2)*ones(K, 1), locU(:, 1) - locUAV(1)*ones(K, 1));

    % h_Los = exp(1j*2*pi*dSU/Lambda)';
    h_Los = exp(1j*2*pi*(locS(:, 1) * (cos(theta).*cos(phi))' + locS(:, 2)* (cos(theta).*sin(phi))') / Lambda);
    
    hBar = sqrt(beta) * (sqrt(kappa/(kappa+1))*h_Los + sqrt(1/(kappa+1))*h_NLoS);
    h = hBar .* repmat(sqrt(sqrt(sum((repmat(locUAV, K, 1) - locU).^2, 2)).^(-alpha))', nIRSrow * nIRScol, 1);
end