function [k_dash,q_k,q_k_dash,rand_num] = proposalDensityk(rand_k,k,...
  rand_num,k_min,k_max)

%[k_dash,q_k,q_k_dash,rand_num] = proposalDensityk(rand_k,k,rand_num,...
% k_min,k_max)
%
% returns the new number of layers (k_dash), proposal for the previous
% number of layers (q_k) and new number of layers (q_k_dash) and a random
% integer determining type of move (rand_num)
%
% -------- INPUT VARIABLES ---------
% rand_k    = indicator vector of type of move  (4 x 1)
% k         = current number of layers
% rand_num  = random integer determining type of move
% k_min     = minimum number of layers
% k_max     = maximum number of layers
%
% -------- OUTPUT VARIABLES ---------
% k_dash    = new number of layers
% q_k       = proposal probability for the previous number of layers in the
% reverse move 
% q_k_dash  = proposal probability for the new number of layers in the
% forward move
% rand_num  = random integer determining type of move; same as input, but
% is updated for edge cases to reflect conservation of probability

if k == k_min
  % only birth, perturbation and no perturbation moves are possible  

  if rand_num <= rand_k(2)
    k_dash   = k+1;
    q_k_dash = 8/10;
    q_k      = 4/10;
    rand_num = 2;
  elseif rand_num > rand_k(2)
    k_dash   = k;
    q_k_dash = 1/10;
    q_k      = 1/10;
  end

elseif k == k_min + 1
  % all moves are possible  
    
  if rand_num <= rand_k(1)
    k_dash   = k+1;
    q_k_dash = 4/10;
    q_k      = 4/10;
  elseif rand_num > rand_k(1) && rand_num <= rand_k(2)
    k_dash   = k-1;
    q_k_dash = 4/10;
    q_k      = 8/10; % q_k changes to account for conservation of probability
  elseif rand_num == rand_k(3)
    k_dash   = k;
    q_k_dash = 1/10;
    q_k      = 1/10;
  else
    k_dash   = k;
    q_k_dash = 1/10;
    q_k      = 1/10;
  end

elseif k > k_min + 1 && k < k_max - 1
  % all moves are possible  

  if rand_num <= rand_k(1)
    k_dash   = k+1;
    q_k_dash = 4/10;
    q_k      = 4/10;
  elseif rand_num > rand_k(1) && rand_num <= rand_k(2)
    k_dash   = k-1;
    q_k_dash = 4/10;
    q_k      = 4/10;
  elseif rand_num == rand_k(3)
    k_dash   = k;
    q_k_dash = 1/10;
    q_k      = 1/10;
  else
    k_dash   = k;
    q_k_dash = 1/10;
    q_k      = 1/10;
  end

elseif k == k_max - 1
  % all moves are possible  

  if rand_num <= rand_k(1)
    k_dash   = k+1;
    q_k_dash = 4/10;
    q_k      = 8/10; % q_k changes to account for conservation of probability
  elseif rand_num > rand_k(1) && rand_num <= rand_k(2)
    k_dash   = k-1;
    q_k_dash = 4/10;
    q_k      = 4/10;
  elseif rand_num == rand_k(3)
    k_dash   = k;
    q_k_dash = 1/10;
    q_k      = 1/10;
  else
    k_dash   = k;
    q_k_dash = 1/10;
    q_k      = 1/10;
  end

elseif k == k_max
  % only death, perturbation and no perturbation moves are possible  
  if rand_num <= rand_k(2)
    k_dash   = k-1;
    q_k_dash = 8/10;
    q_k      = 4/10;
    rand_num = 5;
  elseif rand_num > rand_k(2)
    k_dash   = k;
    q_k_dash = 1/10;
    q_k      = 1/10;
  end

end

end