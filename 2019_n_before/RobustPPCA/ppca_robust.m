function [W, mu, sigma2, C, M, v] = ppca_robust(X, q, v);

%
% Robust probabilistic PCA
%
%  
% args:
%  X -- data
%  q -- number of PCs to retain
%  v -- degree of freedom in t distribution (initial value)
%
% returns:
%  W -- loading matrix
%  mu -- mean
%  sigma2 -- residual variance
%  C -- covariance matrix in data space
%  M -- covariance matrix in score space
%  v -- degree of freedom
%

b_debug = 1;
[N, d] = size(X);

%% initialize W, mu and sigma2

if b_debug % fix random numbers for debugging
    randn('seed', 1); rand('seed', 1);
end;

W = randn(d, q);
mu = mean(X)';
sigma2 = 1;

matId = eye(d); matIq = eye(q);

%% EM algorithm

niter = 100;
L_old = -1e100;
p=zeros(N,1);

for i = 1 : niter;
    
    mat_sig = sigma2*matId;
    C = W * W' + mat_sig;
    isig = 1/sigma2; mat_isig = isig*matId;
    invM = inv ( sigma2*eye(q,q) + W' * W );
    
    [invC, logdetC] = woodbury(mat_sig, W, matIq, W', 1, mat_isig, matIq, 1);
    
    %% calculate log-likelihood up to a constant

    L = -0.5 * N * logdetC; L1 = 0;
    E = X' - repmat(mu, 1, N);
    for j = 1 : N;
        e = E(:,j);
        p(j) = e' * invC * e;
        L1 = L1 + log ( 1 + p(j)/v );
    end
    L = L - 0.5 * L1 * (v+d); L = L / N;
    fprintf('Iter %d, Average LH = %f\n', i, L);

    %    if ( i==10 )    
    if ( L - L_old < 1e-2 )
        break;
    else
        L_old = L;
    end

    
    %% stage 1
    
    % <u_n> & <u_n> x_n
    u = 0; logu = 0;
    u_n = zeros(N,1);
    u_x = zeros(d,1);
    for j = 1 : N;
        u_n(j) = (v+d) / (v+p(j));
        logu_n = psi( 0.5*(v+d) ) - log( 0.5*(v+p(j)) );
        logu = logu + logu_n;
    end
        
    % mu & v
    u_n = u_n / sum(u_n);
    mu = X' * u_n;
    v = fminbnd(@(v) loglik_df(v, N, u, logu), 1, 100);
    
    %% stage 2
    
    E = X' - repmat(mu, 1, N);
    
    if q/d > 0.5;        
        u = 0;
        S = zeros(d,d);

        for j = 1 : N;
            e = E(:,j);
            p(j) = e' * invC * e;
            u_n(j) = (v+d) / (v+p(j));
            u = u + u_n(j);
        end
        S = (E*diag(u_n)*E')/N;
        W1 = S * W * pinv( sigma2*eye(q,q) + invM*W'*S*W );
        sigma2 = trace ( S - S*W*invM*W1' ) / d;
        W = W1;
    else
        u_e = 0;
        t_n = zeros(q,N);        u_t_n = t_n;
        u_tt = zeros(q,q);       u_tt_n = zeros(q,q,N);
        invMW = invM*W';
        t_n = invMW * E;
        sig_invM = sigma2*invM;
        for j = 1 : N;
            e = E(:,j);
            p(j) = e' * invC * e;
            u_n(j) = (v+d) / (v+p(j));
            u = u + u_n(j);
            u_e = u_e + u_n(j)*e'*e;
        end
        u_t_n = repmat(u_n',q,1) .* t_n;
        for j = 1 : N;
            u_tt_n(:,:,j) = sig_invM + u_t_n(:,j)*t_n(:,j)';
        end
        u_tt = sum(u_tt_n,3);
        W = (E*u_t_n') * inv(u_tt);
        
        tr1 = 0; tr2 = 0;
        WE = W'*E; WW = W'*W;
        for j = 1 : N;
            tr1 = tr1 + u_t_n(:,j)' * WE(:,j);
            tr2 = tr2 + trace(WW*u_tt_n(:,:,j));
        end
        sigma2 = (u_e - 2*tr1 + tr2)/(N*d);
    end

end

[R, D] = eig(W'*W);
W = W * R;
C = W * W' + sigma2 * eye(d, d);
M = W' * W + sigma2 * eye(q, q);

