function [u, G, g, h] = updateDeployment(Nt, K, Ns, R, GBar, gBar, hBar, W, thetaMat, tCurrent, muCurrent, locCurrent, locU, locJam, locBS, alpha, P, P2)

    cvx_begin quiet
    cvx_solver mosek;
        variable t 
        variable mu_Opt(K, 1)
        variable omega
        variable omegaBar
        variable varpi(R, 1)
        variable xi(K, 1)
        variable xiBar(K, 1)
        variable o(K, 1)
        variable eOpt(K, R)
        variable u(3, 1)

        maximize t;
        subject to
            [omegaCurrent, omegaBarCurrent, varpiCurrent, Lambda, normASQ, eta, normBSQ] = ComputeParametersSet1(locCurrent, locJam, locBS, alpha, R);
            pow_p(omega, -2/alpha) <= Lambda' * (u - locBS) - normASQ;
            % --- C5f ---
            lhs = norm(u - locBS, 2) + omegaBarCurrent^(-1/alpha - 1) * omegaBar / alpha - (1 + 1/alpha) * omegaBarCurrent^(-1/alpha);
            lhs <= 0;

            pow_p(varpi, -2/alpha) <= diag(eta' * (repmat(u, 1, R) - locJam)) - normBSQ;
            % pow_p(varpiCurrent, -2/alpha)
            % diag(eta' * (repmat(locCurrent, 1, R) - locJam)) - normBSQ
            [xiCurrent, xiBarCurrent, psi, normDSQ] = ComputeParametersSet2(K, locCurrent, locU, alpha);

            % --- RIS Power ---
            lhsPow = 1e6 * P2;
            rhsPow = 0.5*(lhsPow - 1);
            lhsPow = 0.5*(lhsPow + 1);

            rhsPow = [rhsPow; omega * norm(thetaMat*GBar*W, 'fro')/1e2; sqrt(P) * reshape(varpi .* sum(abs(thetaMat * gBar), 1).', R, 1)/1e2; norm(diag(thetaMat), 2)/1e2];
            norm(rhsPow, 2) <= lhsPow;

            % --- Deployment ---
            u(1) >= 0; u(1) <= 100;
            u(2) >= 0; u(2) <= 100;
            u(3) >= 60; u(3) <= 60;

            % 
            % abs(u(1) - locCurrent(1)) <= 0.5;
            % abs(u(2) - locCurrent(2)) <= 0.5;
            muCurrent = muCurrent / 1e5;

            for k = 1 : K
                % --- C5a ---
                lambda = tCurrent / (muCurrent(k));
                lhsA = 1e5 * xiBar(k) * abs(hBar(:, k)'*thetaMat*GBar*W(:, k))^2 + omegaBar;
                rhsA = [2 * sqrt(lambda/2) * mu_Opt(k), 2 * sqrt(1/(2*lambda)) * t, 1e5 * xiBar(k) * abs(hBar(:, k)'*thetaMat*GBar*W(:, k))^2 - omegaBar];
                lhsA >= norm(rhsA, 2);

                % a = 1e5 * xiBarCurrent(k) * abs(hBar(:, k)'*thetaMat*GBar*W(:, k))^2 + omegaBarCurrent;
                % b = [2 * sqrt(lambda/2) * muCurrent(k), 2 * sqrt(1/(2*lambda)) * tCurrent, 1e5 * xiBarCurrent(k) * abs(hBar(:, k)'*thetaMat*GBar*W(:, k))^2 - omegaBarCurrent];
                % norm(b, 2)
                
                % --- C5e ---
                lhsA1 = norm(u - locU(:, k), 2) + xiBarCurrent(k)^(-1/alpha - 1) * xiBar(k) / alpha - (1 + 1/alpha) * xiBarCurrent(k)^(-1/alpha);
                lhsA1 <= 0;

                % a = norm(locCurrent - locU(:, k), 2) + xiBarCurrent(k)^(-1/alpha - 1) * xiBarCurrent(k) / alpha - (1 + 1/alpha) * xiBarCurrent(k)^(-1/alpha);

                % --- C5b ---
                lhsB = 1e10*o(k) * norm(hBar(:, k)'*thetaMat*GBar*W(:, setdiff(1:K, k)), 2)^2 + 1e10*P * sum(eOpt(k, :) .* abs(hBar(:, k)'*thetaMat*gBar).^2) + xi(k)*norm(hBar(:, k)'*thetaMat, 2)^2 + 1;
                rhsB = 1e5*mu_Opt(k);

                lhsB <= rhsB;

                % --- C5c ---
                lambdaC = xiCurrent(k) / omegaCurrent;
                
                rhsC = 1e8*o(k);

                lhsC = 0.5*(rhsC - 1);
                rhsC = 0.5*(rhsC + 1);

                lhsC = [lhsC; sqrt(1e8*lambdaC/2) * (omega); sqrt(1e8*1/(2*lambdaC)) * (xi(k))];

                norm(lhsC, 2) <= rhsC;

                % --- C5d ---
                lambdaD = xiCurrent(k) ./ varpiCurrent;

                rhsD = 1e8*eOpt(k, :).';
                lhsD = 0.5*(rhsD - 1);
                rhsD = 0.5*(rhsD + 1);
                
                lhsD = [lhsD, sqrt(1e8*lambdaD/2) .* (varpi), sqrt(1e8*1./(2*lambdaD)) * (xi(k))];
                
                norms(lhsD, 2, 2) <= rhsD;

                % --- C5g ---
                lhsG = pow_p(xi(k), -2/alpha);
                rhsG = psi(:, k)' * (u - locU(:, k)) - normDSQ(k);

                lhsG <= rhsG;
                % pow_p(xiCurrent(k), -2/alpha)
                % psi(:, k)' * (locCurrent - locU(:, k)) - normDSQ(k)
            end


    cvx_end

    if ~any(isnan(u))
        u = u;
    else
        u = locCurrent;
    end

    G = sqrt(norm(u' - locBS')^(-alpha)) *  GBar;

    % --- Jam to RIS channel
    g = gBar .* sqrt(sqrt(sum((repmat(u', R, 1) - locJam').^2, 2)).^(-alpha))';

    % --- RIS to user channel
    h = hBar .* sqrt(sqrt(sum((repmat(u', K, 1) - locU').^2, 2)).^(-alpha))';

end


function [omegaCurrent, omegaBarCurrent, varpiCurrent, Lambda, normASQ, eta, normBSQ] = ComputeParametersSet1(locCurrent, locJam, locBS, alpha, R)

    omegaCurrent = norm(locCurrent - locBS, 2)^(-alpha);
    omegaBarCurrent = omegaCurrent;
    varpiCurrent = (sqrt(sum((repmat(locCurrent, 1, R) - locJam).^2, 1)).^(-alpha)).';

    Lambda = 2 * (locCurrent - locBS);
    normASQ = sum((locCurrent - locBS).^2);

    eta = 2 * (repmat(locCurrent, 1, R) - locJam);
    normBSQ = (sum((repmat(locCurrent, 1, R) - locJam).^2, 1)).';

end


function [xiCurrent, xiBarCurrent, psi, normDSQ] = ComputeParametersSet2(K, locCurrent, locU, alpha)

    xiCurrent = (sqrt(sum((repmat(locCurrent, 1, K) - locU).^2, 1)).^(-alpha)).';

    xiBarCurrent = xiCurrent;

    psi = 2 * (repmat(locCurrent, 1, K) - locU);
    normDSQ = sum((repmat(locCurrent, 1, K) - locU).^2, 1);

end