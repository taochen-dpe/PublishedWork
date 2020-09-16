%
% calculate matrix inversion & log determinant using woodbury indentity
%
% ( A + BCD )^{-1} = A^{-1} - A^{-1} B ( C^{-1} + D A^{-1} B )^{-1} D A^{-1}
% 
% 
% Ref: 
%  1. Mardia KV, Kent JT & Bibby JM, Multivariate Analysis, Academic
%      Press, 1999
%  2. http://en.wikipedia.org/wiki/Matrix_determinant_lemma
%
function [iV, logdet] = woodbury(A, B, C, D, bdet, iA, iC, bdiagA)

if nargin < 5; bdet = 0; end;

if nargin < 6
    iA = inv(A);
    iC = inv(C);
end

if nargin < 8; bdiagA = 0; end;

Q = iC + D*iA*B;
iV = iA - iA * B * inv(Q) * D * iA;

if bdet
    if bdiagA;
        a = sum(log(diag(A)));
    else
        a = det(A); if a<1e-80; warning('matrix close to singular'); a=1e-80; end;
        a = log(a);
    end
    b = det(Q); if b<1e-80; warning('matrix close to singular'); b=1e-80; end;
    c = det(C); if c<1e-80; warning('matrix close to singular'); c=1e-80; end;
    logdet = a + log(b) + log(c);
end