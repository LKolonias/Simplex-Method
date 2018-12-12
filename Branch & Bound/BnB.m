function [x,cx] = BnB(A,a,c)
%-- Implementation of the Branch & Bound algorithm for 
%   Integer Linear Programming.
%
%   Inputs:
%   A: m x n matrix of inequalities (>=)
%   a: m x 1 right hand side
%   c: optimization function
%
%   Outputs:
%   x: best integer solution
%   cx: cost of the integer solution due to c

[m,n] = size(A);

% Run Simplex algorithm to find the best solution for
% the "relaxed" linear problem
[x,cx,y] = simplex(A,a,c);

if ~isempty(x)
    
    % Create binary tree of the subproblems
    for i=1:length(x)
        if mod(x(i),1)~=0
            e = zeros(1,n);
            e(i)=1;
            
            % Right node of tree (Problem P1)
            P1 = [A;e];
            p1 = [a;ceil(x(i))];
            
            e(i)=-1;
            % Left node of tree (Problem P2)
            P2 = [A;e];
            p2 = [a;-floor(x(i))];
            
            % Call BnB recursively until it returns empty solution
            [x1,cx1] = BnB(P1,p1,c)
            [x2,cx2] = BnB(P2,p2,c)
            
            % Define which of the two solutions has best cost
            if ~isempty(cx1) && ~isempty(cx2)
                if cx1>cx2 
                    x = x1;     cx = cx1;
                else
                    x = x2;     cx = cx2;
                end
            end
            
            if isempty(cx1) && ~isempty(cx2)
                x = x2;     cx = cx2;
            end
            if ~isempty(cx1) && isempty(cx2)
                x = x1;     cx = cx1;
            end
            
            % if there is already integer solution for the specific x_i
            % go to the next one else break
            break;
        else
            
        end
    end

% Return empty solution if the problem is a leaf of the binary tree
else
    x = []; cx = []; 
end