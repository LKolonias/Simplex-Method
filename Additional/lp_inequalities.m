function [h, F, f] = lp_inequalities(n, c, A, a, B, b, J)
%--- lp_inequalities converts a LPG problem in a LP problem
%    with inequalities
% INPUTS
%    n:  number of x's
%    c:  optimization function
%    A:  matrix with equalities
%    a:  right hand side of A
%    B:  matrix with inequalities
%    b:  right hand side of B
%    J:  line vector of 0's and 1's where 1 declares 
%        that the variable is signed
%
% OUTPUTS
%    h:  updated optimization function
%    F:  matrix with all the inequalities
%    f:  right hand side of D

%- Construction of D and d
F = [B;A;-A];
f = [b;a;-a];

%- Put signs in all unsigned variables
F = [F,-F(:,J==0)];
h = [c; -c(J==0)];