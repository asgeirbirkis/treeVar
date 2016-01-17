% Solve the Lorenz equations with treeVar

% Define anonymous functions for differential equations and initial conditions:
odeFun = @(t,u,v,w) [diff(u) - 10*(v - u);
    diff(v) - u.*(28 - w) + v;
    diff(w) - u.*v + (8/3)*w];
icFun = @(u,v,w) [w - 20 ; v + 15; u + 14];

% Right hand side and timespan
rhs = [0; 0; 0];
odeDom = [0 50];

% Solve and plot solution compomentents against each other
[t, uvw] = treeVar.solveIVP(odeFun, icFun, rhs, odeDom);
u = uvw(:,1); v = uvw(:, 2); w = uvw(:, 3);
plot3(u, v, w)