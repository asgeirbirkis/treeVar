% Solve the van der Pol equations with treeVar

% Define anonymous functions for differential equations and initial conditions:
odeFun = @(t,u) diff(u, 2) - 10*(1-u.^2).*diff(u) + u;
icFun = @(u) [u - 1; diff(u)];

% Right hand side and timespand
rhs = 0;
odeDom = [0 100];

% Solve and plot solution
[t, y] = treeVar.solveIVP(odeFun, icFun, rhs, odeDom);
plot(t, y(:, 1))