function pass = test_printTree(~)
%%TEST_PRINTTREE   Do TREEVAR computations, check that print and printTree works.

%% Basic computation:
u = treeVar();
v = cos(u);
w = sin(u);
t = v + w;

%% Call the print method
pass(1) = doesNotCrash(@() print(t));
%% Introducing differentiation
u = treeVar();
myfun = @(u) 2 + diff(u,2);
t = myfun(u);
pass(2) = doesNotCrash(@() print(t));
%% Nested differentiation
t = diff(diff(u)) + diff(u) + u;
pass(3) = doesNotCrash(@() print(t));
end


function pass = doesNotCrash(fn)
try
    fn();
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end
end