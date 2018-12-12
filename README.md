# Simplex Algorithm in MATLAB

Implementation of Simplex algorithm with Lexicographical minimum ratio test and Intiger pivoting for minimization. Solves linear problems of the form LP>= <br><br>

<b>Inputs:</b><br>
  A: m x n matrix of inequalities (>=)<br>
  a: m x 1 right hand side<br>
  c: optimization function<br>
<br>
<b>Outputs:</b><br>
  x: best feasible solution<br>
  cx: cost of the solution due to c<br>
  y: solution of the dual problem


# Additional functions

<b>lp_inequalities:</b> converts a LPG problem in a LP problem with inequalities<br>
<b>lp_equalities:</b> converts a LPG problem in a LP problem with equalities<br>

# Branch & Bound

Implementation of the Branch & Bound algorithm for Integer Linear Programming.

<b>Inputs:</b><br>
  A: m x n matrix of inequalities (>=)<br>
  a: m x 1 right hand side<br>
  c: optimization function<br>
<br>
<b>Outputs:</b><br>
  x: best integer solution
  cx: cost of the integer solution due to c
