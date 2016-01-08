t = treeVar([0 0]);
x = treeVar([1 0]);
y = treeVar([0 1]);
f = diff(x,2) + exp(y) + cos(2*t)
plot(f)