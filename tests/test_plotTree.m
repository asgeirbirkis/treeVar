function pass = test_plotTree(~)
%TEST_PLOTTREE   Do TREEVAR computations, check that PLOTTREE works.

% Create a hidden figure, so that we don't get figures popping up while running
% the test.
hfig = figure('Visible', 'off');

%% Basic computation:
u = treeVar();
v = cos(u);
w = sin(u);
t = v + w;

% Call the plotting method
pass(1) = doesNotCrash(@() plot(t));

%% Introducing differentiation
u = treeVar();
myfun = @(u) 2 + diff(u,2);
s = myfun(u);
pass(2) = doesNotCrash(@() plot(s));

%% Nested differentiation
s2 = diff(diff(u)) + diff(u) + u;
pass(3) = doesNotCrash(@() plot(s2));
end


function pass = doesNotCrash(fn)
try
    fn();
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end
end
