function initVals = setupICs(icFun, diffOrders)
% SETUPICS    Setup initial conditions for ODEs
%
%   Calling sequence:
%       INITVALS = SETUPICS(ICFUN, TOTALDIFFORDERS)
%   where the inputs are
%       ICFUN:      An anonymous function that specifies the initial conditions
%                   of the problem.
%       DIFFORDERS: A vector of the highest differential orders that appear for
%                   each variable in the problem.
%   and the output is
%       INITVALS:   A vector of initial values for variables in their
%                   derivatives in the problem being solved.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check how many unknowns appear in ICFUN.
numArgs = nargin(icFun);
args = cell(numArgs, 1);

% The ID vector to be passed to the TREEVAR constructor.
argsVec = zeros(1, numArgs);

% Populate the args cell with TREEVAR objects.
for argCount = 1:numArgs
    % Set the ID of the current variable to 1:
    argsVec(argCount) = 1;
    % Construct the TREEVAR:
    args{argCount} = treeVar(argsVec);
    % Reset the index vector:
    argsVec = 0*argsVec;
end

% Evaluate ICFUN with the TREEVAR arguments:
icResults = icFun(args{:});

% Initialise vector for initial conditions:
initVals = zeros(1, length(icResults));

for icCounter = 1:length(icResults)
    
    % Current tree we're looking at:
    icTree = icResults(icCounter);
    
    % Go throught the current tree, remove the branch of the tree where the
    % unknown variable appears, and call the methods in the tree with their
    % scalar arguments. The constant that remains is the initial condition
    % specified by the current tree:
    initVals(icCounter) = removeVar(icTree);
        
end

% Obtain indices for initial conditions (i.e. what condition belongs to u and
% what condition belongs to u'):
idx = treeVar.sortConditions(icFun, diffOrders);

% Sort initVals according to the indices
initVals = -initVals(idx);

end

function treeOut = removeVar(treeIn)
% REMOVEVAR    Throw away the leafs where the dependent variable appear

if ( ~isa(treeIn, 'treeVar') )
    % Scalar input, return the same
    treeOut = treeIn;
elseif ( treeIn.height == 0 || ...
        ( treeIn.height == 1 && strcmp(treeIn.method,'diff')))
    % Replace any appearances of unknown variables with 0
    treeOut = 0;
else
    % We're dealing with a + or -
    leftTree = removeVar(treeIn.left);   %#ok<NASGU>
    rightTree = removeVar(treeIn.right); %#ok<NASGU>
    
    % Evaluate the method with arguments that remain. The output will be a
    % scalar that gets propagated back.
    treeOut = eval([treeIn.method,'(leftTree, rightTree)']);
end

end