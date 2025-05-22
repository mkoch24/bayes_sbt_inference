function [p_z_k] = priorProbzk(n,k,r)

% Determine the prior probability of z given k

n_div = n-2*r;     % remove first r and last r intervals
n_rk  = n_div-k*r; % n_rk = n-rk

i = 0;
p_z_k = 1;
while r-i > 0
    p_z_k = p_z_k*(n_rk+r-i)/(n_rk+r-k-i); 
    i = i+1;
end
p_z_k = p_z_k*(1/(n_rk-k));

end