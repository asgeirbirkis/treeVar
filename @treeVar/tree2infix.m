function [infixOut, varArray] = tree2infix(tree, eqno, indexStart)
%TREE2INFIX   Convert a syntax tree to infix form.
%   Calling sequence:
%      [OUT, VARCOUNTER, VARARRAY] = TREE2INFIX(TREE, EQNO, INDEXSTART)
%   where the inputs are
%      TREE:       A syntax tree.
%      EQNO:       An integer, denoting the number of the equation in a coupled
%                  system we're currently converting.
%      INDEXSTART: A vector that denotes at which index we should start indexing
%                  each variable from so that they're in the correct order for
%                  the MATLAB ODE solvers.
%   and the outputs are
%      INFIXOUT:   A string, which describes TREE on infix form.
%      VARARRAY:   A cell, whose first column consist of strings of unique 
%                  identifiers of variables that appear in the equation we're
%                  currently converting to infix form.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If the TREE is empty, do nothing.
if isempty(tree)
    infixOut = '';
    return
end

% Initialize a counter for variables that appear in the current equation.
varCounter = 1;

% Initialize an empty VARARRAY, which will contain identifiers and values of
% variables that appear in the output expression.
varArray = [];

% Start the recursive call for converting a syntax tree to infix form.
[infixOut, varCounter, varArray] = ...
    toInfix(tree, eqno, indexStart, varCounter, varArray);
end


function [out, varCounter, varArray] = ...
    toInfix(tree, eqno, indexStart, varCounter, varArray)
%TOINFIX   Recursively convert syntax tree of infix form.

% Check whether we're converting a syntax tree to infix form, or whether we're
% working directly with a scalar or an anonymous function. If the TREE input
% variable is not a TREEVAR we've arrived at a leaf of the tree, which will be a
% variable. Indicate this by setting NUMARGS = -1, which allows us to use the
% switch statement below.
if ( isa(tree, 'treeVar') )
    % Store the number of arguments.
    numArgs = tree.numArgs;
else
    % Denote that we're at a scalar by setting NUMARGS = -1.
    numArgs = -1;
end

% Convert the current tree to infix form, depending on the number of arguments
% of the operator at the current level.
switch numArgs
    case -1
        % We're looking a a leaf of the syntax tree that is either a variable.
        % Generate a string so that we can include it as a variable later in the
        % conversion process.
        [out, varArray, varCounter] = ...
            scalarInfix(eqno, varCounter, tree, varArray);
    case 0
        % We've hit the TREEVAR constructor level. This might either be the
        % independent variable t, or a dependent variable u (with the
        % appropriate index). Return the appropriate string, based on the ID of
        % the tree:
        if ( ~any(tree.ID) )
            out = 't';
        else
            out = sprintf('u(%i)', indexStart(tree.ID));
        end
    case 1
        % Unary operator tree, convert it recursively to a infix form string.
        [tempOut, varCounter, varArray] = ...
            toInfix(tree.center, eqno, indexStart, varCounter, varArray);
        % Generate a string using the information obtained above, and return it.
        out = sprintf('%s(%s)', tree.method, tempOut);
    case 2
        % Binary operator tree.
        
        if ( strcmp(tree.method, 'diff') )
            % We've hit differentiation. Generate a string, containing the name
            % of the dependent variable -- u --, the differential order, stored
            % in TREE.RIGHT, and index the U string using information from the
            % INDEXSTART vector. Note that TREE.ID is a Boolean vector, so we'll
            % automatically pick out the correct index via INDEXSTART(TREE.ID).
            out = sprintf('u(%i)', tree.right + indexStart(tree.ID));
            return
        end
        
        if ( isa(tree.left, 'treeVar') )
            % Left child node is a tree, convert it recursively to infix form.
            [leftInfix, varCounter, varArray] = ...
                toInfix(tree.left, eqno, indexStart, varCounter, varArray);
        else
            % Left child node is a scalar, generate an infix string we can use
            % later in the process.
            [leftInfix, varArray, varCounter] = ...
                scalarInfix(eqno, varCounter, tree.left, varArray);
        end
        
        if ( isa(tree.right, 'treeVar') )
            % Right child node is a tree, convert it recursively to infix form.
            [rightInfix, varCounter, varArray] = ...
                toInfix(tree.right, eqno, indexStart, varCounter, varArray);
        else
            % Left child node is a scalar, generate an infix string we can use
            % later in the process.
            [rightInfix, varArray, varCounter] = ...
                scalarInfix(eqno, varCounter, tree.right, varArray);
        end
        
        % Combine the infix forms of the left and right nodes into one string on
        % infix form we can return.
        out = sprintf('%s(%s, %s)', tree.method, ...
            leftInfix, rightInfix);
end

end

function [infixOut, varArray, varCounter] = ...
    scalarInfix(eqno, varCounter, tree, varArray)
%SCALARCHEBFUNINFIX   Convert a scalar expression to infix form.
%   A scalar expression is one that only contains a scalar variable.

% The variable name contains which equation we're converting (so that we can
% load the correct variables in the workspace later on), as well as the index
% of variable in the current equation.
varName = sprintf('eq%i_var%i', eqno, varCounter);
infixOut = varName;

% Add the name and value of the current variable to the VARARRAY cell.
varArray = [varArray; {varName, tree}];

% Increase the counter of variables appearing in the equation.
varCounter = varCounter + 1;
end
