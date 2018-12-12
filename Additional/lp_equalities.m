function [g, D, d] = lp_equalities(n, c, A, a, B, b, J)
%--- lp_equalities converts a LPG problem in a LP problem
%    with equalities
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
%    g:  updated optimization function
%    D:  matrix with all the inequalities
%    d:  right hand side of Dble is signed

%- size of matrices
[k, ~] = size(B);
[m,~] = size(A);

%- Constraction of F and f
D = [A, zeros(m,k); B, -eye(k)];
d = [a;b];

%- Update c and J after the slack variables
c = [c; zeros(k,1)];
J = [J, ones(1,k)];

%- Put signs in all unsigned variables
D = [D,-D(:,J==0)];
g = [c; -c(J==0)];