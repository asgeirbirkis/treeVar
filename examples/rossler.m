% Solve the Rossler equations with treeVar

% Define anonymous functions for differential equations and initial conditions:
odeFun = @(t,u,v,w) [diff(u) + v + w; 
    diff(v) - u - 0.2*v; diff(w) - 0.2 - w.*(u - 5.7)];
icFun = @(u, v, w) [u+4; v; w];

% Right hand side and timespan
rhs = [0; 0; 0];
odeDom = [0 100];

% Solve and plot solution compomentents against each other
[t, uvw] = treeVar.solveIVP(odeFun, icFun, rhs, odeDom);
u = uvw(:,1); v = uvw(:, 2); w = uvw(:, 3);
plot3(u, v, w)