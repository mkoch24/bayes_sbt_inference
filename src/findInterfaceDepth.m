function [z,id_z,flag_z] = findInterfaceDepth(z,k,k_dash,...
  rand_num,delta_int,id_int,r,y_depth,rand_k,flag_z)

% z(:,1:2),id_z,flag_z] = find_interface_depth(z(:,1:2),k,k_dash,rand_num,
% delta_int,id_int,r,y_depth,rand_k,flag_z);
% 
% computes the depth of the new interface (birth and perturbation move) and
% returns the new interface depth vector (for all four moves) along with
% the id of the position of the interface in the depth vector which has
% been affected by the move
%
% -------- INPUT VARIABLES --------
% z         = interface depth matrix containing current (initialized to
% zero) and previous sample interface depths (k_max+1 x 2)
% k         = number of layers of previous sample 
% k_dash    = number of layers of current sample 
% rand_num  = random integer determining type of move
% delta_int = interval spacing of the grid in which interfaces can be placed
% id_int    = interval ids where interfaces are present (k+1 x 1)
% r         = minimum interval spacing between two interfaces
% y_depth   = observation data depths (d x 1)
% rand_k    = indicator vector of type of move  (4 x 1)
% flag_z    = flag to ensure that there is always enough space so that a 
% new interface maybe placed (only in the case of birth move); initialized
% to zero
% -------- OUTPUT VARIABLES --------
% z         = interface depth matrix containing current (updated)
% and previous sample interface depths (k_max+1 x 2)
% id_z      = position id of new interface depth in the updated interface
% depth vector z(:,1)
% flag_z    = flag to ensure that there is always enough space so that a 
% new interface maybe placed (only in the case of birth move); initialized
% to zero

z_depth = (0:delta_int:y_depth(end))';

% birth move
if rand_num <= rand_k(1)
  count_z_tmp = 0;
  z_slot_avlbl_tmp = zeros(id_int(end)-1,k);
  for i = 1:k
    if id_int(i+1)-id_int(i) > 2*r+1
      z_slot_avlbl_tmp(1:length(id_int(i)+r+1:id_int(i+1)-r-1),i) = ...
        id_int(i)+r+1:id_int(i+1)-r-1; 
    else
      count_z_tmp = count_z_tmp+1;
    end  
  end

  % check to determine if there is no empty slot to place the new
  % interface. flag_z will always be 0 if k_max is chosen in such a way
  % that there is some space available to place the k_max - 1 interfaces r
  % slots apart
  if count_z_tmp == k
    flag_z  = 1;
    z(:,2)  = z(:,1);
    id_z     = 1;
  end
  
  if flag_z == 0
    z_slot_avlbl = reshape(z_slot_avlbl_tmp,k*(id_int(end)-1),1);
    z_slot_avlbl(z_slot_avlbl == 0) = [];

    z_new_index = z_slot_avlbl(randi([1 length(z_slot_avlbl)],1,1),1);
    z_new       = z_depth(z_new_index)+rand*(z_depth(z_new_index+1)-...
      z_depth(z_new_index));
    z_check     = [z(1:k,1)  z(2:k+1,1)];
    id_z         = find(z_new > z_check(:,1) & z_new < z_check(:,2)); %id of slot in which new interface is placed in z   
    z(id_z+1,2)  = z_new;
    z([1:id_z id_z+2:k_dash+1]',2) = z(1:k+1,1);
  end

% death move
elseif rand_num > rand_k(1) && rand_num <= rand_k(2)

  id_z            = randi([2 k],1); %id of interface to be deleted
  z_temp          = z(1:k+1,1);
  z_temp(id_z)    = [];
  z(1:k,2)        = z_temp;

% perturb interface
elseif rand_num == rand_k(3)

  id_z            = randi([2 k],1);  %id of interface to be perturbed
  z(:,2)          = z(:,1);    
  z(id_z,2)       = z_depth(id_int(id_z),1)+rand*...
    (z_depth(id_int(id_z)+1)-z_depth(id_int(id_z)));

% no change of interface
else

  id_z            = [];
  z(1:k_dash+1,2) = z(1:k+1,1);

end

end