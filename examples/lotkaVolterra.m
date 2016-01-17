% Solve the Lotka-Volterra predator prey equations with treeVar

% Define anonymous functions for differential equations and initial conditions:
odeFun = @(t,u, v) [diff(u) - 2*u + u.*v; diff(v) + v - u.*v];
icFun = @(u,v) [u-1; v-.5];

% Right hand side of DEs and timespan
rhs = [0; 0];
odeDom = [0 100];

% Call the solver!
[t, uv] = treeVar.solveIVP(odeFun, icFun, rhs, odeDom);

% Plot the solutions
plot(t, [uv(:, 1) uv(:, 2)])