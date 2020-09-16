%
% Given mean (mu) and covariance matrix (C), work out the 
%  distribution of x as Gaussian if there's missing data present
%

function [m, S] = prob_miss(x, mu, C)

m = x;  S = zeros(size(C));
idmis = find(isnan(x));
idobv = find(~isnan(x));

if (isempty(idmis))
    return;
end

m1 = mu(idobv);
m2 = mu(idmis);
x1 = x(idobv);

C11 = C(idobv, idobv);
C12 = C(idobv, idmis);
C21 = C(idmis, idobv);
C22 = C(idmis, idmis);

W = C21 * pinv(C11);
x2 = m2 + W * (x1 - m1);
cov2 = C22 - W * C12;

m(idmis) = x2;
S(idmis, idmis) = cov2;