function [WOpt, thetaMat] = updateBeamMinMax(Nt, K, Ns, R, thetaMatCurrent, wCurrent, dCurrent, tCurrent, muCurrent, G, g, h, Pt, P2, P, alpha)
    
    % ------ channels
    d = (h' * thetaMatCurrent * G)';
    cvx_begin quiet
        cvx_solver mosek;
        variable W(Nt, K) complex;
        variable t;

        maximize t;
        subject to
            % --- Transmit Power ---
            sum_square_abs(vec(W)) <= Pt;

            % --- RIS Power ---
            sum_square_abs(vec(thetaMatCurrent*G*W)/1e2) + P * norm(thetaMatCurrent*g, 'fro')^2/1e4 + norm(diag(thetaMatCurrent), 2)^2/1e4 <= 1e6*P2;
            
            for k = 1 : K
                [bCurrent, ACurrent] = ComputeParametersSet1(k, dCurrent, wCurrent, thetaMatCurrent, g, h, P);

                lhsA = bCurrent - real(ACurrent(1, 1)) - t - real(ACurrent(2, 2)) * (P * norm(h(:, k)' * thetaMatCurrent * g, 2)^2 + norm(h(:, k)' * thetaMatCurrent, 2)^2 + 1);
                rhsA = 2*real(ACurrent(1, 2) * dCurrent(:, k)' * W(:, k)) + real(ACurrent(2, 2)) * sum_square_abs(vec(dCurrent(:, k)' * W));

                lhsA >= rhsA;
            end

    cvx_end

    if ~any(isnan(W(:)))
        WOpt = W;
    else
        % WOpt = wCurrent;
        WOpt = d * (d' * d)^(-1);
        WOpt = sqrt(Pt) * WOpt / norm(WOpt);
    end

    cvx_begin quiet
        cvx_solver mosek;
        variable theta(Ns, 1) complex;
        variable tOpt;

        maximize tOpt;
        subject to
            abs(theta) <= alpha; % amp

            % --- RIS Power ---
            sum_square_abs(vec(diag(theta)*G*WOpt)/1e2) + P * sum_square_abs(vec(diag(theta)*g)/1e2) + sum_square_abs(theta/1e2) <= 1e6*P2;
            
            for k = 1 : K
                [bCurrent, ACurrent] = ComputeParametersSet1(k, dCurrent, wCurrent, thetaMatCurrent, g, h, P);

                lhsA = bCurrent - real(ACurrent(1, 1)) - tOpt - real(ACurrent(2, 2));
                rhsA = 2*real(ACurrent(1, 2) * (h(:, k)' * diag(theta) * G) * WOpt(:, k)) + real(ACurrent(2, 2)) * sum_square_abs(vec((h(:, k)' * diag(theta) * G) * WOpt)) + real(ACurrent(2, 2)) * (P * sum_square_abs(vec(h(:, k)' * diag(theta) * g)) + sum_square_abs(vec(h(:, k)' * diag(theta))));

                lhsA >= rhsA;
            end

    cvx_end

    if ~any(isnan(theta))
        thetaMat = diag(theta);
    else
        thetaMat = thetaMatCurrent;
    end

end


function [bCurrent, ACurrent] = ComputeParametersSet1(k, dCurrent, wCurrent, thetaMatCurrent, g, h, P)

    D = [1; 0];
    
    BCurrent = [1, wCurrent(:, k)' * dCurrent(:, k);
                dCurrent(:, k)' * wCurrent(:, k), dCurrent(:, k)' * (wCurrent * wCurrent') * dCurrent(:, k) + P * norm(h(:, k)' * thetaMatCurrent * g, 2)^2 + norm(h(:, k)' * thetaMatCurrent, 2)^2 + 1];

    ACurrent = BCurrent^(-1) * D * (D' * BCurrent^(-1) * D)^(-1) * D' * BCurrent^(-1);

    bCurrent = real(log(D' * BCurrent^(-1) * D) + trace(ACurrent * BCurrent));

end