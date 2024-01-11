function [xhat1, xhat2, xhat3, xhat4] = l2ARsmoothing2(y, lam, a, b)
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
    N = length(y)+2;
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
        Upsilon = Upsilon + diag(a(i) * ones(N - abs(i - 1), 1), i - 1);
    end
    Upsilon = Upsilon(1:N-3,:);
    Upsilon = sparse(Upsilon);

    Gamma = zeros(N);
    for i = 1:length(b)
        Gamma = Gamma + diag(b(i) * ones(N - abs(i - 1), 1), i - 1);
    end
    Gamma = Gamma(1:N-2,:);
    Gamma = sparse(Gamma);
    % GG = inv(Gamma*Gamma');
    % Gright = Gamma'*GG;
    % Gleft = GG*Gamma;
    % D = Upsilon*Gright;
    % [U,S,V] = svd(D);
    % A = eye(N-2) + lam*S'*S;
    % 
    % M = Gright*V*inv(A)*V'*Gleft;
    % xhat1 = M*(Gamma' * Gamma*y);
    % X = chol(inv(A),'lower');
    % L = X*GG*Gamma;
    %xhat2 = Gamma'*V*L*y;

    % GG = inv(Gamma*Gamma');
    % D = Upsilon*Gamma'*GG;
    % A = speye(N-2) + lam*D'*D;
    % xhat1 = A\y;
    % L = chol(A\eye(N-2),'lower');
    % xhat2 = L'*y; %xhat2 = L*xhat2;

    M = Gamma' * Gamma + lam * Upsilon'* Upsilon;
    xhat1 = Gamma * inv(M)* Gamma' * y;
    L = chol(inv(M),'lower');
    xhat2 = L*Gamma'* y;
    xhat2 = xhat2(2:end-1);
    % GG = inv(Gamma*Gamma');
    % F = Gamma'*GG*GG*Gamma*Upsilon'*inv(1/lam*eye(N-3) + Upsilon*Gamma'*GG*GG*Gamma*Upsilon')*Upsilon;
    % xhat2 = y - F*y;
    % L = chol(inv(1/lam*eye(N-3) + Upsilon*Gamma'*GG*GG*Gamma*Upsilon'), 'lower');
    % xhat2 = y - L*y;

    % D = Upsilon*Gamma'*GG;
    % A = speye(N-2) + lam*D'*D;
    % xhat1 = A\y(2:end-1);
    % L = chol(A\eye(N-2),'lower');
    % xhat2 = L*y(2:end-1);


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