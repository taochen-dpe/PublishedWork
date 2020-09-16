%
% 
%

function value = loglik_df (v, N, e_u, e_logu)

value = 0.5*N*v*log(v/2) - N*gammaln(v/2) + (0.5*v-1)*e_logu - 0.5*v*e_u;
value = - value;