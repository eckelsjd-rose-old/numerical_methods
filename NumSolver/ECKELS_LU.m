
% The function YOURLASTNAME_LU takes in a 3x3 matrix, A, and outputs the LU
% decomposition matrices, L and U. This function does NOT perform partial
% pivoting using permutation matrices.
% 
% Inputs:
%           A - 3x3 matrix
% Outputs:
%           L - 3x3 lower triangular matrix
%           U - 3x3 upper triangular matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Change the function name AND the file name with your last name.  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,U] = ECKELS_LU(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%        Write the code to compute the matrices L and U.        %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l21 = A(2,1) / A(1,1);
l31 = A(3,1) / A(1,1);
A(2,:) = A(2,:) - l21*A(1,:);
A(3,:) = A(3,:) - l31*A(1,:);
l32 = A(3,2) / A(2,2);
A(3,:) = A(3,:) - l32*A(2,:);
U = A;
L = [1 0 0; l21 1 0; l31 l32 1];

