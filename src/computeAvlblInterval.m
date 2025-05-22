function [n_avlbl] = computeAvlblInterval(id_int,r)

% n_avlbl = compute_avlbl_interval(id_int,r)
%
% determine the number of intervals available to place new interface such
% that new interface is at least r intervals apart from existing interfaces
%
% -------- INPUT VARIABLES ---------
% id_int    = interval ids where interfaces are present (k+1 x 1)
% r         = min interval separation
% -------- OUTPUT VARIABLES --------
% n_avlbl   = number of intervals available to place new interface

id_int_diff        = diff(id_int);
id_int_diff(id_int_diff <= 2*r) = [];

if isempty(id_int_diff)
  n_avlbl  = 0;
else  
  n_avlbl   = sum(id_int_diff-2*ones(length(id_int_diff),1).*r-1);    
end

end