function [W, mu, sigma2, C, M, v] = ppca_robust_miss(X, q, v);

%
% Robust probabilistic PCA with missing data
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



[N, d] = size(X);

%% initialize W, mu and sigma2
%randn('seed', 1); rand('seed', 1);
W = randn(d, q);
%W = zeros(d,q);
mu = zeros(d,1);
for i = 1 : d;
    ti = X(:,i);
    mu(i) = mean (ti(find(~isnan(ti))));
end
sigma2 = 1;

%% EM algorithm

niter = 100;
L_old = -1e100;

for i = 1 : niter;
    
    C = W * W' + sigma2 * eye(d, d);
    isig = 1/sigma2;
    invM = pinv ( sigma2*eye(q,q) + W' * W );
    %invC = pinv(C);
    invC = isig*eye(d,d) - W*invM*W'*isig; % Woodbury inversion identity

    %% calculate log-likelihood up to a constant
    L = 0;
    for j = 1 : N;
        x = X(j,:)';
        id_o = find(~isnan(x));
        d_o = length(id_o);
        
        e = x(id_o) - mu(id_o);
        p = e' * invC(id_o,id_o) * e;
        L = L + gammaln((v+d_o)/2) - 0.5 * log(det(C(id_o,id_o))) ...
            - 0.5*d_o*log(pi*v) - gammaln(v/2) - 0.5 * (v+d_o)* log ( 1 + p/v );
    end
    L = L / N;
    fprintf('Iter %d, Average LH = %.2f, df=%.2f\n', i, L, v);
    if ( L - L_old < 1e-1 )
        break;
    else
        L_old = L;
    end

    
    %% stage 1
    
    % <u_n> & <u_n> x_n
    u = 0; logu = 0;
    u_x = zeros(d,1);
    for j = 1 : N;

        x = X(j,:)';
        [z, Q] = prob_miss(x, mu, C);
        
        id_o = find(~isnan(x));
        d_o = length(id_o);
        
        e = x(id_o) - mu(id_o);
        p_o = e' * invC(id_o,id_o) * e;
        u_n = (v+d_o) / (v+p_o);
        logu_n = psi( 0.5*(v+d_o) ) - log( 0.5*(v+p_o) );
        logu = logu + logu_n;        
        u = u + u_n;
        u_x = u_x + u_n*z;
    end
        
    % mu
    mu = u_x / u;
    v = fminbnd(@(v) loglik_df(v, N, u, logu), 1, 100);
    
    %% stage 2
    
    % S
    u = 0;
    S = zeros(d,d);
    for j = 1 : N;
        
        x = X(j,:)';
        [z, Q] = prob_miss(x, mu, C);
        
        id_o = find(~isnan(x));
        d_o = length(id_o);
        
        e = x(id_o) - mu(id_o);        
        p_o = e' * invC(id_o,id_o) * e;
        u_n = (v+d_o) / (v+p_o);
        u = u + u_n;

        e = z - mu;
        S = S + Q + u_n * e * e';
    end
    S = S / N;
    
    W1 = S * W * pinv( sigma2*eye(q,q) + invM*W'*S*W );
    sigma2 = trace ( S - S*W*invM*W1' ) / d;
    W = W1;

end

[R, D] = eig(W'*W);
W = W * R;
C = W * W' + sigma2 * eye(d, d);
M = W' * W + sigma2 * eye(q, q);

