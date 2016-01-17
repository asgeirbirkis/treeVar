% Solve the IVP
%   x''+ Ax + ax^3 + cy^2 = 0
%   y''+ Ay + ay^3 + cx^2 = 0
% with initial conditions:
%   x(0) = 0.7, x'(0) = 0, y(0) = 0.8, y'(0) = 0.
A = 0.3; a = 1; c = 1;  % Specify coefficients
% Define the second order coupled system of ODEs
odeFun = @(t,x,y) [diff(x,2) + A*x + a*x.^3 + c*y.^2; ...
    diff(y,2) + A*y + a*y.^3 + c*x.^2];
% Anonymous function for initial conditions:
icFun = @(x,y) [x - .7; diff(x); y - .8; diff(y)];
% Interval we want to work on
odeDom = [0 100];
% RHS of ODEs
rhs = [0; 0];

% Call the solver!
[t, xy] = treeVar.solveIVP(odeFun, icFun, rhs, odeDom);
% Plot the solutions
plot(t, [xy(:, 1) xy(:,3)])
