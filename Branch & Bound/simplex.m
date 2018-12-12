function [x,cx,y,T] = simplex(A,a,c)
%-- Implementation of Simplex algorithm with Lexicographical
%   minimum ratio test and Intiger pivoting for minimization
%   linear problems of the form LP>=
%
%   Inputs:
%   A: m x n matrix of inequalities (>=)
%   a: m x 1 right hand side
%   c: optimization function
%
%   Outputs:
%   x: best feasible solution
%   cx: cost of the solution due to c
%   y: solution of the dual problem


% ============ PHASE I ====================================================

% Construction of tableau T
[m,n] = size(A);
T = [A,-eye(m);c',zeros(1,m)];

I(1:m) = n+1:(n+m);     % I holds the columns of the base

J = find(-a<0);
for i=1:length(J)
    e = zeros(m+1,1);
    e(J(i))=1;
    T = [T,e];
end

T = [T,[-a;0];zeros(1,n+m),ones(1,length(J)),0;];

% Pivot in of the new variables
for i=1:length(J)
    T = jx(T,J(i),m+n+i);
    I(J(i)) = m+n+i;    
end

[m,n] = size(T);

% Start integer pivoting of phase I with LexMin Ratio
PrevPivot = -1;
while min(T(m,1:(n-1)))<0
    % Find column of pivot element
    s = T(m,1:(n-1));
    [el, col] = min(s);
    
    % Find row of pivot element
    [~,srt] = LexMinRatio(T,col,1);
    CurPivot = T(srt(1),col);
    
    % Store column on base
    I(srt(1)) = col;
    
    % Integer pivoting
    for i=1:m
        if i~=srt(1)
            T(i,:) = (-T(srt(1),col)*T(i,:) + T(i,col)*T(srt(1),:)) / abs(PrevPivot);
        end
    end
   PrevPivot = CurPivot;
   
end

% In the end of phase I check if solution is feasible
% else return empty
if max(T(m,:))>0
    x = []; cx = []; y = [];
    return;
end

% Remove z1 row from tableau
T = T(1:m-1,:);

% ============ PHASE II ===================================================
[m,n] = size(T);

% Run simplex algorithm with Integer pivoting and LexMin Ratio
while min(T(m,1:(n-1)))<0

    s = T(m,1:(n-1));
    s(s>0)=min(s)-1; s(s==0)=min(s)-1;
    [el, col] = max(s);
    [~,srt] = LexMinRatio(T,col,2);
    
    % Check if Linear Problem is not bounded
    if isempty(srt)
        x = []; cx = []; y = [];
        return;   
    end
    
    CurPivot = T(srt(1),col);
    I(srt(1)) = col;
    for i=1:m
        if i~=srt(1)
            T(i,:) = (-T(srt(1),col)*T(i,:) + T(i,col)*T(srt(1),:)) / abs(PrevPivot);
        end
    end
   PrevPivot = CurPivot; 
end

I = sort(I);
d = -T(find(T(:,I(1))),I(1))

% Return minimum cost solution vector
[m2, n2] = size(A);
for i =1:n2
    x(i,1) = T(find(T(:,I(i))),n)/d;
end

% Return cost of the solution
cx = -T(m,n)/d;

% Return the solution of dual problem
y = T(m,setdiff(1:m2+n2,I))'/d;