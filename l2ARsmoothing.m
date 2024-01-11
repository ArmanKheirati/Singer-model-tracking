function [xhat1, xhat2, xhat3, xhat4] = l2ARsmoothing(y, lam, a, b)
%%%  ARMA smoothing filter using blockwise matrix formuation and forward-backward 
%%%  filtering 
%%%  y is the noisy input signal
%%%  lam is the regularization factor
%%%  a is the vector containing AR coefficents 
%%%  b is the vector containing MA coefficents
%%%  xhat1 is the smoothed output using block-wise matrix formulation
%%%  xhat2 is the filtered output using FIR block-wise matrix formulation
%%%  xhat3 is the smoothed output using forward filtering
%%%  xhat4 is the smoothed output using forward bacward filtering
    N = length(y);
    b = b/norm(b);
    a = a/norm(a);
    if a(1) == 0 
        disp('Error: $\beta = 0$ means it is not a singer model')
        a = a(2:end);
    end
    if b(1) == 0
        disp('Error: $\beta = 1$ means it is not a singer model')
        b = b(2:end);
    end
    if ~isreal(b)
        disp('Error: $\zeta$ is not real')
        b = abs(b);
    end

    Upsilon = zeros(N);
    for i = 1:length(a)
        Upsilon = Upsilon + diag(a(i) * ones(N - i + 1, 1), i - 1);
    end
    Upsilon = Upsilon(1:N-3,:);
    Upsilon = sparse(Upsilon);

    Gamma = zeros(N);
    for i = 1:length(b)
        Gamma = Gamma + diag(b(i) * ones(N - i + 1, 1), i - 1);
    end
    Gamma = Gamma(1:N-2,:);
    Gamma = sparse(Gamma);

    M = Gamma' * Gamma + lam * Upsilon'* Upsilon;
    xhat1 = inv(M)*(Gamma' * Gamma*y);
    L = chol(inv(M),'lower');
    xhat2 = L*Gamma'* y(2:size(Gamma,1)+1);
    % %%%%%%%%%%%%% Woudbary matrix Identity  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % G1 = inv(Gamma*Gamma');
    % GG = Gamma'*G1*G1*Gamma;
    % 
    % UUT = Upsilon * GG * Upsilon';
    % Uy = Upsilon*y;
    % F = sparse(1:size(Upsilon,1), 1:size(Upsilon,1), 1/lam) + UUT; % F : Sparse banded matrix
    % xhat1 = y - GG*Upsilon'*(F\Uy); % Solve banded linear system
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%  Forward-backward filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phi = conv(a, a(end:-1:1));
    gg = phi*lam;

    psi = conv(b, b(end:-1:1));
    psi = [0 psi 0];
    gg = gg + psi;

    r = roots(gg);
    r_abs = abs(r);
    I_causal = r_abs < 1;
    r_causal = r(I_causal);
    p_causal = poly(r_causal);

    psi = conv(b, b(end:-1:1));
    r = roots(psi);
    r_abs = abs(r);
    I_causal = r_abs < 1;
    r_causal = r(I_causal);
    b_causal = poly(r_causal);
    K = sum(p_causal)/sum(b_causal);
    xhat3 = filtfilt(K*b_causal, p_causal, y);
    xhat4 = filter(K*b_causal, p_causal, y);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end