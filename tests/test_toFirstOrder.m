function pass = test_toFirstOrder()
%TEST_TOFIRSTORDER   Test conversion to first-order format for various DEs.

%% Setup
tol = 5e-14;

%% Simple test
% We have the expression 5*(diff(u, 2) + 3*u), so expect the first order
% system to be:
%   u'(1) = u(2)
%   u'(2) = -3*u(1)
problemNo = 1;
myfun = @(x, u) 5*(diff(u, 2) + 3*u);
rhs = 0;
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
pass(1, problemNo) = norm(anonFun(1, [2 1]) - [1;-6]) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = norm(coeffs{1} - 5) < tol;
pass(4, problemNo) = all( diffOrders == 2);

%% Try with variables in the anonymous function:
problemNo = 2;
alpha = 4;
myfun = @(x, u) 3.5*(diff(u, 2) + alpha*u);
rhs = 5;

% Expect the first order system to be:
%   u'(1) = u(2)
%   u'(2) = rhs/3.5 - alpha*u(1)
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
pass(1, problemNo) = norm(anonFun(1,[2 1]) - [1; 5/3.5-4*2]) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = norm(coeffs{1} - 3.5) < tol;
pass(4, problemNo) = all( diffOrders == 2);

%% Try with a variable coefficient in the anonymous function:
problemNo = 3;
myfun = @(x, u) 5*((x+1).*diff(u, 2) + u);
rhs = -5;
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
x = -.5;
% Expect this to be [u(2); (-5-5*u(1))/(5*(x(-.5)+1))] = [1; -6]
pass(1, problemNo) = norm( anonFun(-.5,[2 1]) - [1;-6] ) < tol;
pass(2, problemNo) = (idx == 1 );
pass(3, problemNo) = norm(feval(coeffs{1}, x) - 5*(x+1)) < tol;
pass(4, problemNo) = all( diffOrders == 2);

%% Another expression with a variable coefficient in the expression
problemNo = 4;
myfun = @(x,u) diff(u, 2) + 2*sin(x).*u;
x = .5;
rhs = cos(x);
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);

% Expect this to be [u(2); cos(x(.5)) - 2*sin(x(.5))]
pass(1, problemNo) = norm( anonFun(.5,[2 1]) - [1; cos(.5)-2*(sin(.5)*2)]) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = ( norm(coeffs{1} - 1) < tol);
pass(4, problemNo) = all( diffOrders == 2);

%% Variable coefficient at the start
problemNo = 5;
x = -.5;
myfun = @(x,u) 3*(x+2).*((x+1).*diff(u, 2) + u);
rhs = 0;
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);

% Expect this to be [u(2); -u(1))/(x(-.5)+1)] = [1; -4]
pass(1, problemNo) = norm( anonFun(-.5,[2 2]) - [2; -4]) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = ( norm(feval(coeffs{1},x) - 3*(x+2).*(x+1)) < tol);
pass(4, problemNo) = all( diffOrders == 2);

%% Piecewise smooth problem
problemNo = 6;
x = .5;
myfun = @(x,u) cos(x).*diff(u, 2) + 2*abs(sin(pi*x)).*u;
rhs = cos(2*x);
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) [u(2); (cos(2*x)-2*abs(sin(pi*x)).*u(1))./cos(x)];
pass(1, problemNo) = norm(anonFun(.5,[2 1]) - correctFun(.5, [2 1])) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = norm(feval(coeffs{1}, x) - cos(x)) < tol;
pass(4, problemNo) = all( diffOrders == 2);

%% Piecewise smooth, higher order derivatives
problemNo = 7;
t = .5;
myfun = @(x,u) cos(x).*diff(u, 4) + 2*abs(sin(pi*x)).*diff(u,2);
rhs = cos(2*x);
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) [u(2); u(3); u(4); ...
    (cos(2*x)-2*abs(sin(pi*x)).*u(3))./cos(x)];
pass(1, problemNo) = norm(anonFun(.5,[2 1 3 4]) - correctFun(.5, [2 1 3 4])) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = norm(feval(coeffs{1},t) - cos(t)) < tol;
pass(4, problemNo) = all( diffOrders == 4);

%% Simple coupled system, first order
% We have the equations 
%   5*(diff(u) + 3*v) = tanh(x)
%   cos(x)*diff(v) + sin(x).*u = 2
% so expect the first order system to be:
%   u'(1) = tanh(x)/5 - 3*u(2)
%   u'(2) = (2 - sin(x)*u(1))/cos(x)
problemNo = 8;
t = 1;
x = 1;
myfun = @(x,u,v) [5*(diff(u) + 3*v); cos(x).*diff(v) + sin(x).*u];
rhs = [tanh(x); 2];
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) [tanh(x)/5 - 3*u(2); (2-sin(x).*u(1))/cos(x)];
evalPt = [2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == [1 2]);
pass(3, problemNo) = norm(coeffs{1} - 5) + ...
    norm(feval(coeffs{2},t) - cos(x)) < tol;
pass(4, problemNo) = all( diffOrders == [1 1]);

%% Simple coupled system, second order
% We have the equations 
%   5*(diff(u,2) + 3*v) = tanh(x)
%   cos(x)*diff(v,2) + sin(x).*u = 2
% so expect the first order system to be:
%   u'(1) = u(2)
%   u'(2) = tanh(x)/5 - 3*u(3)
%   u'(3) = u(4)
%   u'(4) = (2 - sin(x)*u(1))/cos(x)
problemNo = 9;
myfun = @(x,u,v) [5*(diff(u,2) + 3*v); cos(x).*diff(v,2) + sin(x).*u];
rhs = [tanh(x); 2];
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) [u(2); tanh(x)/5 - 3*u(3); u(4); (2-sin(x).*u(1))/cos(x)];
evalPt = [2 1 2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == [1 3]);
pass(3, problemNo) = norm(coeffs{1} - 5) + norm(feval(coeffs{2}, t) - cos(x)) < tol;
pass(4, problemNo) = all( diffOrders == 2);

%% Coupled system, second order, piecewise smooth
% We have the equations 
%   5*(diff(u,2) + 3*v) = tanh(x)
%   cos(x)*diff(v,2) + abs(sin(pi*x)).*u = 2
% so expect the first order system to be:
%   u'(1) = u(2)
%   u'(2) = tanh(x)/5 - 3*u(3)
%   u'(3) = u(4)
%   u'(4) = (2 - sin(x)*u(1))/cos(x)
problemNo = 10;
myfun = @(x,u,v) [5*(diff(u,2) + 3*v); 
    cos(x).*diff(v,2) + abs(sin(pi*x)).*u];
rhs = [tanh(x); 2];
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) [u(2); tanh(x)/5 - 3*u(3); u(4); (2-abs(sin(pi*x)).*u(1))/cos(x)];
evalPt = [2 1 2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == [1 3]);
pass(3, problemNo) = norm(coeffs{1} - 5) + norm(feval(coeffs{2}, t) - cos(x)) < tol;
pass(4, problemNo) = all( diffOrders == 2);

%% Variable coefficient inside the anonymous function, rather than on the RHS:
% We have the equations 
%   5*(diff(u,2) + coth(x+5)) = -5
% so expect the first order system to be:
%   u'(1) = u(2)
%   u'(2) = -5/5 - coth(x+5)
problemNo = 11;
myfun = @(x,u) 5*(diff(u, 2) + coth(x+5));
rhs = -5;
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) [u(2); -1 - coth(x+5)];
evalPt = [2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = (idx == 1 );
pass(3, problemNo) = norm(coeffs{1} - 5) < tol;
pass(4, problemNo) = all( diffOrders == 2);

%% Four variables, mixed derivatives, mixed order
% We have the equations
%   diff(y,3) + w + diff(u) = exp(x)
%   7*diff(w) + diff(v) = x.^2 
%   5*(diff(u,2) + 3*v + coth(x)) = tanh(x)
%   cos(x)*diff(v,2) + diff(u) - diff(y,2) = 2
% so expect the first order system to include the variables:
%   u(1) = u, u(2) = u', u(3) = v, u(4) = v'; u(5) = w, u(6) = y, u(7) = y',
%   u(8) = y''
% and the first order reformulation to be
%   u'(1) = u(2)
%   u'(2) = tanh(x)/5 - 3*u(3) - coth(x+5)
%   u'(3) = u(4)
%   u'(4) = (2 - u(2) + u(8))/cos(x)
%   u'(5) = (x.^2-u(4))/7
%   u'(6) = u(7)
%   u'(7) = u(8)
%   u'(8) = exp(x) - u(5) - u(2)
problemNo = 12;
myfun = @(x,u,v,w,y) [...
    diff(y, 3 ) + w + diff(u);
    7*diff(w) + diff(v);
    5*(diff(u,2) + 3*v + coth(x+5));
    cos(x).*diff(v,2) + diff(u) - diff(y,2)];
rhs = [exp(x); x.^2; tanh(x); 2];
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) [...
    u(2); tanh(x)/5 - 3*u(3) - coth(x+5); 
    u(4); (2-u(2) + u(8))/cos(x);
    (x.^2 - u(4))/7;
    u(7); u(8); exp(x) - u(5) - u(2)];
evalPt = [2 1 2.4 2.3 1.2 3.2 5.1 4.6];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == [1 3 5 6]);
pass(3, problemNo) = norm(coeffs{1}(t) - 5) + norm(coeffs{2}(t) - cos(x)) + ...
    norm(coeffs{3}(t) - 7) + norm(coeffs{4}(t) - 1) < tol;
pass(4, problemNo) = all( diffOrders == [2 2 1 3]);

%% Scalar problem, multiple variable coefficients and scalars in operator, ver1
% We have the equations 
%   diff(u) + x + 3*x = tanh(x)
% so expect the first order system to be:
%   u'(1) = tanh(x) - 4*x
problemNo = 13;
myfun = @(x,u) diff(u) + x + 3*x;
rhs = tanh(x);
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) tanh(x) - 4*x;
evalPt = [2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == 1);
pass(3, problemNo) = norm(coeffs{1}(t) - 1) < tol;
pass(4, problemNo) = all( diffOrders == 1);

%% Scalar problem, multiple variable coefficients and scalars in operator, ver2
% We have the equations 
%   cos(x)*(2 - sin(x) + diff(u)) + x + 3*x = tanh(x)
% so expect the first order system to be:
%   u'(1) = (tanh(x) - 4*x)/cos(x) + sin(x) - 2
problemNo = 14;
myfun = @(x,u) cos(x).*(2 - sin(x) + diff(u)) + x + 3*x;
rhs = tanh(x);
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) (tanh(x) - 4*x)/cos(x) + sin(x) - 2;
evalPt = [2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == 1);
pass(3, problemNo) = norm(coeffs{1}(t) - cos(x)) < tol;
pass(4, problemNo) = all( diffOrders == 1);

%% Scalar problem, multiple variable coefficients and scalars in operator, ver3
% We have the equations 
%   2 - sin(x) + diff(u) + x + 3*x = tanh(x)
% so expect the first order system to be:
%   u'(1) = tanh(x) - 4*x + sin(x) - 2
problemNo = 15;
myfun = @(x,u) 2 - sin(x) + diff(u) + x + 3*x;
rhs = tanh(x);
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) tanh(x) - 4*x + sin(x) - 2;
evalPt = [2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == 1);
pass(3, problemNo) = norm(coeffs{1}(t) - 1) < tol;
pass(4, problemNo) = all( diffOrders == 1);

%% Simple coupled system, multiple variable coefficients and scalars in operator
% We have the equations 
%   4*(2 + diff(u)) + x + 3*x + 3*v = tanh(x)
%   cos(x)*diff(v) + sin(x).*u + x + x = 2
% so expect the first order system to be:
%   u'(1) = (tanh(x) - 4*x - 3*u(2))/4 - 2
%   u'(2) = (2 - 2*x -sin(x)*u(1))/cos(x)
problemNo = 16;
myfun = @(x,u,v) [...
    4*(2 + diff(u)) + x + 3*x + 3*v; 
    cos(x).*diff(v) + sin(x).*u + x + x];
rhs = [tanh(x); 2];
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) [...
    (tanh(x) - 4*x - 3*u(2))/4 - 2;
    (2 - 2*x - sin(x)*u(1))/cos(x)];
evalPt = [2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == [1 2]);
pass(3, problemNo) = norm(coeffs{1}(t) - 4) + norm(coeffs{2}(t) - cos(x)) < tol;
pass(4, problemNo) = all( diffOrders == [1 1]);

%% Introduce breakpoints in the RHS
problemNo = 17;
x = .5;
t = x;
myfun = @(x,u) cos(x).*diff(u, 2) + 2*abs(sin(pi*x)).*u;
rhs = 0;
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) [u(2); ((abs(x/2-round(x/2))<.05)-2*abs(sin(pi*x)).*u(1))./cos(x)];
pass(1, problemNo) = norm(anonFun(.5,[2 1]) - correctFun(.5, [2 1])) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = norm(coeffs{1}(t) - cos(x)) < tol;
pass(4, problemNo) = all( diffOrders == 2);

%% Coupled system, second order, breakpoints in the RHS as well
% We have the equations 
%   5*(diff(u,2) + 3*v) = tanh(x)
%   cos(x)*diff(v,2) + abs(sin(pi*x)).*u = (abs(x/2-round(x/2))<.05)
% so expect the first order system to be:
%   u'(1) = u(2)
%   u'(2) = tanh(x)/5 - 3*u(3)
%   u'(3) = u(4)
%   u'(4) = ((abs(x/2-round(x/2))<.05) - sin(x)*u(1))/cos(x)
problemNo = 18;
x = 1;
t = x;
myfun = @(x,u,v) [5*(diff(u,2) + 3*v); 
    cos(x).*diff(v,2) + abs(sin(pi*x)).*u];
rhs = [tanh(x); 0];
[anonFun, idx, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs);
correctFun = @(x,u) [u(2); tanh(x)/5 - 3*u(3); u(4); ...
    ((abs(x/2-round(x/2))<.05) -abs(sin(pi*x)).*u(1))/cos(x)];
evalPt = [2 1 2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == [1 3]);
pass(3, problemNo) = norm(coeffs{1}(t) - 5) + norm(coeffs{2}(t) - cos(x)) < tol;
pass(4, problemNo) = all( diffOrders == 2);

%% Coupled systems -- Unsupported format, highest order derivatives in same eqn
myfun = @(x,u,v) [diff(u,2) + diff(v,2); diff(u) + sin(v)];
rhs = [1;2];
try
    treeVar.toFirstOrder(myfun, rhs);
    errorPass(1) = 0;
catch ME
    % The highest order derivatives of u and v appear in the same line -- this
    % should give us an error.
    errorPass(1) = strcmp(ME.identifier, 'TREEVAR:toFirstOrder:diffOrders');
end

%% Coupled systems -- Nonlinearity in highest order derivative
myfun = @(x,u,v) [diff(u,2).*diff(v); diff(v,2) + sin(u)];
rhs = [1;2];
try
    treeVar.toFirstOrder(myfun, rhs);
    errorPass(2) = 0;
catch ME
    % We're multiplying the highest order derivative by a variable, this should
    % give an error.
    errorPass(2) = strcmp(ME.identifier, 'TREEVAR:expandTree:nonlinearity');
end

%% Combine the information
pass = [pass(:)' errorPass];
end
