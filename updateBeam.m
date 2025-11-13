function [WOpt, thetaMat, t] = updateBeam(Nt, K, Ns, R, thetaMatCurrent, wCurrent, dCurrent, tCurrent, muCurrent, G, g, h, Pt, P2, P, alpha)
    

    cvx_begin quiet
    cvx_solver mosek;
        variable W(Nt, K) complex
        variable theta(Ns, 1) complex
        variable t
        variable mu_Opt(K, 1)
        variable betaOpt(K, K)
        variable betaOptBar(K, K)
        variable varthetaOpt(K, R)
        variable varthetaOptBar(K, 1)
        variable varkappaOpt(K, Ns)
        variable varkappaOptBar(K, Ns)

        maximize t

        % ------ channels
        d = (h' * diag(theta) * G)';
        C = diag(theta) * G;
        % gN = h' * diag(theta) * g; % K*R

        subject to
            sum_square_abs(vec(W)) <= Pt; % power
            
            % ------ RIS Power ------
            lhsRISPow = 1e6 * P2;
            rhsRISPow = 0.5 * (lhsRISPow - 1);
            lhsRISPow = 0.5 * (lhsRISPow + 1);
            rhsRISPow = [rhsRISPow; reshape(varkappaOpt, K*Ns, 1)/1e2; reshape(varkappaOptBar, K*Ns, 1)/1e2; sqrt(P) * reshape(repmat(theta, 1, R) .* g, Ns*R, 1)/1e2; theta/1e2];

            norm(rhsRISPow, 2) <= lhsRISPow;

            abs(theta) <= alpha; % amp

            [uCurrent, uAbsSQ, uVecCurrent, uCurrentNormSQ, delta] = ComputeParametersSet1(dCurrent, wCurrent);

            for k = 1 : K
                % --- C5a
                lhsA = real(delta(:, k)' * d(:, k) + uVecCurrent(:, k)' * W(:, k)) ... 
                       - 0.5*uCurrentNormSQ(k) - uAbsSQ(k);
                
                rhsA = 0.5 * (lhsA - 1);
                lhsA = 0.5 * (lhsA + 1);
                
                lambda = tCurrent / (muCurrent(k));
                % a = real(delta(:, k)' * dCurrent(:, k) + uVecCurrent(:, k)' * wCurrent(:, k)) ... 
                %        - 0.5*uCurrentNormSQ(k) - uAbsSQ(k);
                % 
                % b = 0.5*(muCurrent(k)^2 * lambda + tCurrent^2/lambda) + 0.5*norm(uCurrent(k)*dCurrent(:, k) - wCurrent(:, k), 2)^2
                rhsA = [rhsA; sqrt(lambda/2) * mu_Opt(k); sqrt(1/(2*lambda)) * t; sqrt(1/2) * (uCurrent(k)*d(:, k) - W(:, k))];
                lhsA >= norm(rhsA, 2);
                
                % --- C5b
                lhsB = 0.5 * mu_Opt(k);

                rhsB = [betaOpt(k, setdiff(1:K, k)), betaOptBar(k, setdiff(1:K, k)), varthetaOpt(k, :), varthetaOptBar(k), lhsB - 1];

                lhsB >= norm(rhsB, 2);

                [Lambda, normCSQ, eta, normDSQ, psi, normESQ, phi, normFSQ] = ComputeParametersSet2(k, K, dCurrent, wCurrent);

                % --- C5c ---
                % imag(diag(repmat(dCurrent(:, k), 1, K-1)'*wCurrent(:, setdiff(1:K, k))))
                % % % 
                % 0.25*(norm(dCurrent(:, k) - 1j*wCurrent(:, 3), 2)^2 + normESQ(2) - 2*real(psi(:, 2)' * (dCurrent(:, k) + 1j*wCurrent(:, 3))))
                % 0.25*(norm(dCurrent(:, k) + 1j*wCurrent(:, 3), 2)^2 + normFSQ(2) - 2*real(phi(:, 2)' * (dCurrent(:, k) - 1j*wCurrent(:, 3))))

                lhsC = reshape(betaOpt(k, setdiff(1:K, k)), K - 1, 1) + 0.5 * real(diag(Lambda' * (repmat(d(:, k), 1, K-1) - W(:, setdiff(1:K, k)))));

                rhsC = 0.5 * (lhsC - 1 - 0.25 * normCSQ);
                lhsC = 0.5 * (lhsC + 1 - 0.25 * normCSQ);

                rhsC = [rhsC, 0.5 * (repmat(d(:, k), 1, K-1) + W(:, setdiff(1:K, k))).'];

                lhsC >= norms(rhsC, 2, 2);

                lhsD = reshape(betaOpt(k, setdiff(1:K, k)), K - 1, 1) + 0.5 * real(diag(eta' * (repmat(d(:, k), 1, K-1) + W(:, setdiff(1:K, k)))));

                rhsD = 0.5 * (lhsD - 1 - 0.25 * normDSQ);
                lhsD = 0.5 * (lhsD + 1 - 0.25 * normDSQ);

                rhsD = [rhsD, 0.5 * (repmat(d(:, k), 1, K-1) - W(:, setdiff(1:K, k))).'];

                lhsD >= norms(rhsD, 2, 2);

                % --- C5d ---
                lhsE = reshape(betaOptBar(k, setdiff(1:K, k)), K - 1, 1) + 0.5 * real(diag(psi' * (repmat(d(:, k), 1, K-1) + 1j*W(:, setdiff(1:K, k)))));

                rhsE = 0.5 * (lhsE - 1 - 0.25 * normESQ);
                lhsE = 0.5 * (lhsE + 1 - 0.25 * normESQ);

                rhsE = [rhsE, 0.5 * (repmat(d(:, k), 1, K-1) - 1j*W(:, setdiff(1:K, k))).'];

                lhsE >= norms(rhsE, 2, 2);

                lhsF = reshape(betaOptBar(k, setdiff(1:K, k)), K - 1, 1) + 0.5 * real(diag(phi' * (repmat(d(:, k), 1, K-1) - 1j*W(:, setdiff(1:K, k)))));

                rhsF = 0.5 * (lhsF - 1 - 0.25 * normFSQ);
                lhsF = 0.5 * (lhsF + 1 - 0.25 * normFSQ);

                rhsF = [rhsF, 0.5 * (repmat(d(:, k), 1, K-1) + 1j*W(:, setdiff(1:K, k))).'];

                lhsF >= norms(rhsF, 2, 2);
                
                % --- C5e ---
                varthetaOpt(k, :).' >= sqrt(P)*abs(h(:,k)'*diag(theta)*g).';

                % --- C5f ---
                varthetaOptBar(k) >= norm(h(:,k)'*diag(theta), 2);

                % --- RIS Power ---
                [Omega, normGSQ, Sigma, normHSQ, Psi, normISQ, Xi, normJSQ] = ComputeParametersSet3(k, Ns, thetaMatCurrent, G, wCurrent);
                % imag(thetaMatCurrent(2,2)*G(2,:)*wCurrent(:, k))
                % 0.25*(norm((thetaMatCurrent(2,2)*G(2,:))' + 1j*wCurrent(:, k),2)^2 + norm((thetaMatCurrent(2,2)*G(2,:))' - 1j*wCurrent(:, k),2)^2 - 2*real(Xi(:, 2)'*((thetaMatCurrent(2,2)*G(2,:))' - 1j*wCurrent(:, k))))
                
                lhsG = reshape(varkappaOpt(k, :), Ns, 1) + 0.5 * real(diag(Omega' * (C' - repmat(W(:, k), 1, Ns))));

                rhsG = 0.5 * (lhsG - 1 - 0.25 * normGSQ);
                lhsG = 0.5 * (lhsG + 1 - 0.25 * normGSQ);

                rhsG = [rhsG, 0.5 * (C' + repmat(W(:, k), 1, Ns)).'];

                lhsG >= norms(rhsG, 2, 2);

                lhsH = reshape(varkappaOpt(k, :), Ns, 1) + 0.5 * real(diag(Sigma' * (C' + repmat(W(:, k), 1, Ns))));

                rhsH = 0.5 * (lhsH - 1 - 0.25 * normHSQ);
                lhsH = 0.5 * (lhsH + 1 - 0.25 * normHSQ);

                rhsH = [rhsH, 0.5 * (C' - repmat(W(:, k), 1, Ns)).'];

                lhsH >= norms(rhsH, 2, 2);

                lhsI = reshape(varkappaOptBar(k, :), Ns, 1) + 0.5 * real(diag(Psi' * (C' + 1j*repmat(W(:, k), 1, Ns))));

                rhsI = 0.5 * (lhsI - 1 - 0.25 * normISQ);
                lhsI = 0.5 * (lhsI + 1 - 0.25 * normISQ);

                rhsI = [rhsI, 0.5 * (C' - 1j*repmat(W(:, k), 1, Ns)).'];

                lhsI >= norms(rhsI, 2, 2);

                lhsJ = reshape(varkappaOptBar(k, :), Ns, 1) + 0.5 * real(diag(Xi' * (C' - 1j*repmat(W(:, k), 1, Ns))));

                rhsJ = 0.5 * (lhsJ - 1 - 0.25 * normJSQ);
                lhsJ = 0.5 * (lhsJ + 1 - 0.25 * normJSQ);

                rhsJ = [rhsJ, 0.5 * (C' + 1j*repmat(W(:, k), 1, Ns)).'];

                lhsJ >= norms(rhsJ, 2, 2);
            end
            

    cvx_end
    % norm(diag(theta)*G*W, 'fro')^2 + P*norm(diag(theta)*g, 'fro')^2 + norm(theta, 'fro')^2
    if ~any(isnan(theta(:))) && ~any(isnan(W(:)))
        thetaMat = diag(theta);
        WOpt = W;
    else
        thetaMat = thetaMatCurrent;
        % WOpt = wCurrent;
        % ------ channels
        d = (h' * thetaMat * G)';
        WOpt = d * (d' * d)^(-1);
        WOpt = sqrt(Pt) * WOpt / norm(WOpt);
    end

end


function [uCurrent, uAbsSQ, uVecCurrent, uCurrentNormSQ, delta] = ComputeParametersSet1(dCurrent, wCurrent)
    
    uCurrent = diag(dCurrent' * wCurrent);
    uVecCurrent = dCurrent .* uCurrent.' + wCurrent;
    uAbsSQ = abs(uCurrent).^2;    
    uCurrentNormSQ = sum(abs(uVecCurrent).^2, 1);
    delta = uVecCurrent .* uCurrent';

end


function [Lambda, normCSQ, eta, normDSQ, psi, normESQ, phi, normFSQ] = ComputeParametersSet2(k, K, dCurrent, wCurrent)

    % --- constant values for C5c ---
    Lambda = repmat(dCurrent(:, k), 1, K-1) - wCurrent(:, setdiff(1:K, k));
    normCSQ = sum(abs(Lambda).^2, 1)';

    eta = repmat(dCurrent(:, k), 1, K-1) + wCurrent(:, setdiff(1:K, k));
    normDSQ = sum(abs(eta).^2, 1)';

    % --- constant values for C5d ---
    psi = repmat(dCurrent(:, k), 1, K-1) + 1j*wCurrent(:, setdiff(1:K, k));
    normESQ = sum(abs(repmat(dCurrent(:, k), 1, K-1) + 1j*wCurrent(:, setdiff(1:K, k))).^2, 1)';

    phi = repmat(dCurrent(:, k), 1, K-1) - 1j*wCurrent(:, setdiff(1:K, k));
    normFSQ = sum(abs(repmat(dCurrent(:, k), 1, K-1) -1j*wCurrent(:, setdiff(1:K, k))).^2, 1)';

end

function [Omega, normGSQ, Sigma, normHSQ, Psi, normISQ, Xi, normJSQ] = ComputeParametersSet3(k, Ns, thetaMatCurrent, G, wCurrent)
    
    % --- theta ---
    casChannel = thetaMatCurrent * G; % M * Ns

    % --- constant values for RIS Power ---
    Omega = casChannel' - repmat(wCurrent(:, k), 1, Ns); % M * Ns
    normGSQ = sum(abs(Omega).^2, 1)';

    Sigma = casChannel' + repmat(wCurrent(:, k), 1, Ns);
    normHSQ = sum(abs(Sigma).^2, 1)';

    % --- constant values for RIS Power ---
    Psi = casChannel' + 1j*repmat(wCurrent(:, k), 1, Ns);
    normISQ = sum(abs(Psi).^2, 1)';

    Xi = casChannel' - 1j*repmat(wCurrent(:, k), 1, Ns);
    normJSQ = sum(abs(Xi).^2, 1)';

end