function [WOpt, thetaMat] = updateBeamSDR(Nt, K, Ns, R, thetaMatCurrent, wCurrent, dCurrent, tCurrent, muCurrent, G, g, h, Pt, P2, P, alpha)
    

    % cvx_begin quiet
    % cvx_solver mosek;
    %     variable W(Nt, K) complex
    %     variable t
    %     variable mu_Opt(K, 1)
    % 
    %     maximize t;
    % 
    %     % ------ channels
    %     d = (h' * thetaMatCurrent * G)';
    %     subject to
    %         % ------ Transmit Power ------
    %         sum_square_abs(vec(W)) <= Pt; % power
    % 
    %         % ------ RIS Power ------
    %         sum_square_abs(vec(thetaMatCurrent*G*W)/1e2) + P * norm(thetaMatCurrent*g, 'fro')^2/1e4 + norm(diag(thetaMatCurrent) ,2)^2/1e4 <= 1e6*P2;
    % 
    %         [uCurrent, uAbsSQ] = ComputeParametersSet1(dCurrent, wCurrent);
    % 
    %         for k = 1 : K
    %             % --- C5a
    %             lhsA = 2 * real(uCurrent(k)' * d(:, k)' * W(:, k)) - uAbsSQ(k);
    % 
    %             rhsA = 0.5 * (lhsA - 1);
    %             lhsA = 0.5 * (lhsA + 1);
    % 
    %             lambda = tCurrent / (muCurrent(k));
    % 
    %             rhsA = [rhsA; sqrt(lambda/2) * mu_Opt(k); sqrt(1/(2*lambda)) * t];
    %             lhsA >= norm(rhsA, 2);
    % 
    %             % --- C5b
    %             lhsB = 0.5 * mu_Opt(k);
    % 
    %             rhsB = [lhsB - 1, d(:, k)' * W(:, setdiff(1:K, k)), sqrt(P) * h(:, k)' * thetaMatCurrent * g, norm(h(:, k)' * thetaMatCurrent, 2)];
    % 
    %             lhsB >= norm(rhsB, 2);
    %         end
    % 
    % 
    % cvx_end
    % 
    % if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    %     WOpt = W;
    % else
    %     WOpt = wCurrent;
    % end

    % ------ channels
    d = (h' * thetaMatCurrent * G)';
    cvx_begin quiet
    cvx_solver mosek;
        variable R(Nt, Nt) hermitian semidefinite
        variable R1(Nt, Nt, K) complex
        
        variable t
        variable mu_Opt(K, 1)
        maximize t;

        subject to
            for k = 1:K
                R1(:,:,k) == hermitian_semidefinite(Nt);
            end
            
            % ------ Transmit Power ------
            trace(R) <= Pt; % power
            
            % ------ RIS Power ------
            real(trace(G'*(thetaMatCurrent' * thetaMatCurrent)*G * R) / 1e4) + P * norm(thetaMatCurrent*g, 'fro')^2/1e4 + norm(diag(thetaMatCurrent) ,2)^2/1e4 <= 1e6*P2;

            for k = 1 : K
                % --- C5a
                lhsA = real(trace(d(:, k) * d(:, k)' * R1(:,:,k)));
                
                rhsA = 0.5 * (lhsA - 1);
                lhsA = 0.5 * (lhsA + 1);
                
                lambda = tCurrent / (muCurrent(k));
                
                rhsA = [rhsA; sqrt(lambda/2) * mu_Opt(k); sqrt(1/(2*lambda)) * t];
                lhsA >= norm(rhsA, 2);
                
                % --- C5b
                lhsB = mu_Opt(k);

                rhsB = real(trace(d(:, k) * d(:, k)' * (R - R1(:,:,k)))) + P * norm(h(:, k)' * thetaMatCurrent * g, 2)^2 + norm(h(:, k)' * thetaMatCurrent, 2)^2 + 1;

                lhsB >= rhsB;
            end

            R == sum(R1, 3);
            
    cvx_end
    WOpt = zeros(Nt, K);
    if ~any(isnan(R1(:)))
        for k = 1 : K
            WOpt(:, k) = R1(:, :, k) * d(:, k) / sqrt(real(d(:, k)' * R1(:, :, k) * d(:, k)));
        end
    else
        WOpt = d * (d' * d)^(-1);
        WOpt = sqrt(Pt) * WOpt / norm(WOpt);
    end

    [V, D] = eig(conj(diag(thetaMatCurrent))*diag(thetaMatCurrent).');
    [~, maxIndex] = max(diag(D));
    maxEigVec = V(:, maxIndex);

    cvx_begin quiet
    cvx_solver mosek;
        variable thetaMatOpt(Ns, Ns) hermitian semidefinite
        variable tOpt;
        variable mu_Opt(K, 1)

        maximize tOpt - 1e-3*(real(trace(thetaMatOpt)) - real(maxEigVec' * thetaMatOpt * maxEigVec)); % 

        diag(thetaMatOpt) <= alpha^2;

        % --- RIS Power ---
        thetaHatMat = diag(diag(thetaMatOpt));
        real(trace(G * (WOpt * WOpt') * G' * thetaHatMat)) / 1e4 + P * real(trace(g * g' * thetaHatMat)) / 1e4 + real(trace(thetaHatMat)) / 1e4 <= 1e6 * P;

        for k = 1 : K

            lhsA = real(trace(diag(h(:, k)') * G * (WOpt(:, k) * WOpt(:, k)') * G' * diag(h(:, k)) * thetaMatOpt));
                
            rhsA = 0.5 * (lhsA - 1);
            lhsA = 0.5 * (lhsA + 1);
            
            lambda = tCurrent / (muCurrent(k));
            
            rhsA = [rhsA; sqrt(lambda/2) * mu_Opt(k); sqrt(1/(2*lambda)) * tOpt];
            lhsA >= norm(rhsA, 2);

            lhsB = mu_Opt(k);

            rhsB = real(trace(diag(h(:, k)') * G * (WOpt(:, setdiff(1:K, k)) * WOpt(:, setdiff(1:K, k))') * G' * diag(h(:, k)) * thetaMatOpt) + P * trace(diag(h(:, k)') * (g * g') * diag(h(:, k)) * thetaMatOpt) + trace(diag(h(:, k)') * diag(h(:, k)) * thetaMatOpt) + 1);

            lhsB >= rhsB;

        end
    cvx_end

    if ~any(isnan(thetaMatOpt(:)))
        [V, D] = eig(thetaMatOpt);
        [maxEigVal, maxIndex] = max(diag(D));
        maxEigVec = sqrt(maxEigVal) * V(:, maxIndex);

        thetaMat =  diag(maxEigVec');

    else
        thetaMat = thetaMatCurrent;
    end

end


function [uCurrent, uAbsSQ] = ComputeParametersSet1(dCurrent, wCurrent)
    
    uCurrent = diag(dCurrent' * wCurrent);
    uAbsSQ = abs(uCurrent).^2;    

end


function thetaMat = GaussianRandomization(K, Ns, wCurrent, thetaMatOpt, G, g, h, P, P2, alpha)

    L = 1000; % number of Gaussian randomizations

    maxValue = 0;
    [U, Sigma] = eig(thetaMatOpt);
    for l = 1 : L
        r = sqrt(2) / 2 * (randn(Ns, 1) + 1j * randn(Ns, 1));
        v = U * Sigma^(0.5) * r;
        v = v(1 : Ns);
        thetaMatCurrent = diag(v');

        flag = CheckFesible(thetaMatCurrent, wCurrent, G, g, P, P2, alpha);

        dCurrent = (h' * thetaMatCurrent * G)';
        gCurrent = h' * thetaMatCurrent * g;
        [minSINR, ~] = K_update_SINR(wCurrent, dCurrent, gCurrent, thetaMatCurrent, h, P, K);
        if minSINR > maxValue && flag
            thetaMat = thetaMatCurrent;
            maxValue = minSINR;
        end
    end

end


function flag = CheckFesible(thetaMatCurrent, wCurrent, G, g, P, P2, alpha)

    flag = 1;

    if any(abs(thetaMatCurrent) > alpha)
        flag = 0;
    end

    if norm(vec(thetaMatCurrent*G*wCurrent)/1e2, 2)^2 + P * norm(vec(thetaMatCurrent*g)/1e2, 2)^2 + norm(diag(thetaMatCurrent)/1e2, 2)^2 > 1e6*P2
        flag = 0;
    end

end
