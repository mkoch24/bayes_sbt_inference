function [state, stats] = runMCMC(data, const, state, stats)

% [struct(state), struct(stats)] = runMCMC(struct(data), struct(const), 
% struct(state), struct(stats))
% 
% runs all three blocks for n_samp iterations
%
% -------- INPUT VARIABLES --------
% structs: data, const, state, stats. Please see initializeState.m
%
% -------- OUTPUT VARIABLES --------
% updated structs state, stats. Please see initializeState.m

for ip = 2:const.n_samp
  % current sample
  state.ip          = ip;
  % determine type of move: birth/death/perturbation/fixed
  state.rand_num    = randi([1 const.rand_k(end)],1);
  % block 1
  [state, stats]    = sampleFirstBlock(data, const, state, stats);
  % block 2
  [state, coeff]    = sampleSecondBlock(data, const, state);
  % block 3
  [state]           = sampleThirdBlock(const, state, coeff);
end

end