function [L,srt] = LexMinRatio(T,col,f)
%-- Function LexMinRatio sorts the rows of a Tableau T for
%   the Simplex algorithm.
%
%   Inputs:
%   T: Tableau of Simplex algorithm
%   col: Column of pivoting
%   f: Defines if the Ratio Test is for Phase I or Phase II
%   (In first case Tableau consists of extra row)
%
%   Outputs:
%   L: Matrix of sorted rows
%   srt: The sorting vector

[m,n] = size(T);
if(f==1); T(m-1,:)=[]; end
J = zeros(1,(m-1));

% Construction of matrix L with the rows of the tableau which
% participate in the Ratio Test
for i=1:m-1
    if (T(i,col)>=0) J(i)=1; end
    L(i,1) = (1/T(i,col))*(-T(i,n));
    k=2;
    for j=n-m+1:n-1          
            L(i,k)=(1/T(i,col))*(T(i,j));
            k=k+1;
    end
end
L=L';   
[m,n] = size(L);
srt(1:n) = 1:n;
L=L(:,J==0);
srt = srt(J==0);

% BubbleSort Lexicographic Ratio on the matrix L
[m,n] = size(L);
dist = n-1;

for j=1:n-1
    for i=1:dist
        v = L(:,i)-L(:,i+1);
        if v(find(v,1))>0
            ex = L(:,i);         ex2 = srt(i);       
            L(:,i) = L(:,i+1);   srt(i) = srt(i+1); 
            L(:,i+1) = ex;       srt(i+1) = ex2;     
        end
    end
dist = dist-1;
end