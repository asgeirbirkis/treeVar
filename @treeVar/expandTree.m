function newTree = expandTree(tree, maxOrder)
%EXPANDTREE   Convert expressions like 5*(diff(u) + u) to 5*diff(u) + 5*u.
%   The role of the EXPANDTREE method is to split up syntax trees, so that
%   expressions involving the highest derivative appearing in the equation stand
%   alone, that is, not inside parenthesis. For example, EXPANDTREE will convert
%   expressions like 5*(diff(u) + u) to 5*diff(u) + 5*u. This is necessary in
%   order to be able to put the correct expression on the right-hand side when
%   calling MATLAB's built-in solvers.
%
%   Calling sequence:
%      NEWTREE = EXPANDTREE(TREE, MAXORDER)
%   where the inputs are:
%      TREE:       The syntax tree we want to split.
%      MAXORDER:   A vector containing the maximum differential order of each
%                  variable that appear in the problem under consideration.
%   and the output is:
%      NEWTREE:    A syntax tree where the highest order deriative has been 
%                  taken out of any parenthesis where other variables with lower
%                  differential order appear.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(tree, 'treeVar') || tree.height <= 1 || ...
        all(tree.diffOrder < maxOrder) || ...
        ( ~tree.hasTerms && sum(tree.ID) <= 1) || ~any(tree.ID) )
    % We don't need to worry about expanding the tree any of the following
    % conditions apply: 
    %   * The input is not a tree
    %   * It is a very short tree
    %   * Its diffOrder is less than the maxOrder we consider
    %   * It doesn't any operators that split it up in terms (that is, + or -),
    %     and it doesn't include more than one variable (which would imply
    %     nonlinearity, such as u.*w).
    %   * No dependent variables appear in it.
    newTree = tree;
    
elseif ( tree.numArgs == 1 )
    % If we have a unary operator, we expand the center node recursively.
    newTree = expandTree(tree.center, maxOrder);
    
elseif ( any(strcmp(tree.method, {'plus', 'minus'})) )
    % We're at a + or a -, so expand the syntax tree recursively:
    tree.left = treeVar.expandTree(tree.left, maxOrder);
    tree.right = treeVar.expandTree(tree.right, maxOrder);
    newTree = tree;
    
else
    % Must be at .*, ./ or .^.
    
    % If both left and right arguments to the operator are syntax trees where
    % dependent variables appear (i.e. any other variables apart from the
    % independent time variable), and the highest order derivative appear in
    % any, we know that we must have nonlinearities in the highest order
    % derivative, e.g. in the case of a first order IVP, u.*diff(u). This is not
    % supported, so throw an error:
    if ( isa(tree.left, 'treeVar') && isa(tree.right, 'treeVar') && ...
            ( any(tree.diffOrder == maxOrder) && ...
            any(tree.left.ID) && any(tree.right.ID)  ) )
        error('TREEVAR:expandTree:nonlinearity', ...
            ['Nonlinearity in highest order derivative detected.\n' ...
            'Unable to convert to first order format.']);
    end
    
    % Find out what kind of argument we're dealing with for the left argument
    % of the operator
    if ( ~isa(tree.left, 'treeVar') )
        % Left argument is not a syntax tree, must be either a CHEBFUN or a
        % scalar. In either case, the left tree has zero children.
        leftNumArgs = 0;
    else
        % Left tree is indeed a syntax tree, store its number of children.
        leftNumArgs = tree.left.numArgs;
    end
    
    % Find out what kind of argument we're dealing with for the right argument
    % of the operator:
    if ( ~isa(tree.right, 'treeVar') )
        % Right argument is not a syntax tree, must be either a CHEBFUN or a
        % scalar. In either case, the right tree has zero children.
        rightDiffOrder = 0;
        rightNumArgs = 0;
    else
        % Right tree is indeed a syntax tree, store its number of children.
        rightNumArgs = tree.right.numArgs;
    end
    
    % Split up the left tree depending on how many children it has. In case we
    % don't need to split at the current level, we store the differential
    % order and the height of the left tree in the variables LEFTDIFFORDER
    % and LEFTHEIGHT. Finally, we introduce the Boolean variables SPLITLEFT
    % which indicates whether the left tree had to be split:
    if ( leftNumArgs == 0 )
        % We're either dealing with a scalar left tree
        leftTree = tree.left;
        leftDiffOrder = 0*tree.diffOrder;
        leftHeight = 0;
        splitLeft = false;      
    elseif ( ~any(tree.left.ID) )
        % Only independent variables on the left
        leftTree = tree.left;
        leftDiffOrder = 0*tree.diffOrder;
        leftHeight = leftTree.height;
        splitLeft = false;
    elseif ( strcmp(tree.left.method,'diff') && tree.left.height <= 1 )
        % We're differentiating one of the basis variables directly (as that's
        % the only case that tree.right.height <= 1). In this case, no splitting
        % is needed.
        leftTree = tree.left;
        leftDiffOrder = leftTree.diffOrder;
        leftHeight = leftTree.height;
        splitLeft = false;
    elseif ( leftNumArgs == 1 )
        % Here we're dealing with the unary operator case. In this case, we
        % don't need to split at the current level, however, we recursively
        % split up the child tree.
        leftTree = treeVar.expandTree(tree.left, maxOrder);
        leftDiffOrder = leftTree.diffOrder;
        leftHeight = leftTree.height;
        splitLeft = false;
    else
        % Here we're dealing with the binary operator case. In this case, we do
        % need to split at the current level, as well as recursively splitting
        % up both the child trees. LEFTLEFT and LEFTRIGHT will be the syntax
        % trees stored in the left child of the current level, once we've split
        % up the child trees of the current level.
        leftLeft = treeVar.expandTree(tree.left.left, maxOrder);
        leftRight = treeVar.expandTree(tree.left.right, maxOrder);
        
        % Indicate that the left tree had to be split
        splitLeft = true;
    end
    
    % Split up the right tree depending on how many children it has. In case we
    % don't need to split at the current level, we store the differential order
    % and the height of the right tree in the variables RIGHTDIFFORDER and
    % RIGHTHEIGHT. Finally, we introduce the Boolean variables SPLITRIGHT which
    % indicates whether the right tree had to be split:
    if ( rightNumArgs == 0 )
        % We're either dealing with a scalar right tree
        rightTree = tree.right;
        rightDiffOrder = 0*tree.diffOrder;
        rightHeight = 0;
        splitRight = false;
    elseif ( ~any(tree.right.ID) )
        % Only independent variables on the right
        rightTree = tree.right;
        rightDiffOrder = 0*tree.diffOrder;
        rightHeight = rightTree.height;
        splitRight = false;
    elseif ( strcmp(tree.right.method,'diff') && tree.right.height <= 1 )
        % We're differentiating one of the basis variables directly (as that's
        % the only case that tree.right.height <= 1). In this case, no splitting
        % is needed.
        rightTree = tree.right;
        rightDiffOrder = rightTree.diffOrder;
        rightHeight = rightTree.height;
        splitRight = false;
    elseif ( rightNumArgs == 1 )
        % He're we're dealing with the unary operator case. In this case, we
        % don't need to split at the current level, however, we recursively
        % split up the child tree.
        rightTree = treeVar.expandTree(tree.right, maxOrder);
        rightDiffOrder = rightTree.diffOrder;
        rightHeight = rightTree.height;
        splitRight = false;
    else
        % He're we're dealing with the binary operator case. In this case, we do
        % need to split at the current level, as well as recursively splitting
        % up both the child trees. RIGHTLEFT and RIGHTRIGHT will be the syntax
        % trees stored in the right child of the current level, once we've split
        % up the child trees of the current level.
        rightLeft  = treeVar.expandTree(tree.right.left, maxOrder);
        rightRight = treeVar.expandTree(tree.right.right, maxOrder);
        
        % Indicate that the right tree had to be split
        splitRight = true;
    end
    
    % If we had to split either the left tree or the right tree, construct new
    % syntax trees which we'll then glue together below. Note that we should
    % never expect to split both the left and right trees, as that indicates
    % nonlinearities in the highest order derivative.
    if ( ~splitLeft && ~splitRight )
        % Neither tree had to be split, so we can simply return the input tree!
        newTree = tree;
        return
        
    elseif ( splitLeft && ~splitRight )
        % Had to split on the left, not the right.
        
        % Multiply together the new left factor of the left tree, and the
        % current right tree:
        newLeft = treeVar();
        newLeft.method = 'times';
        newLeft.numArgs = 2;
        newLeft.left = leftLeft;
        newLeft.right = rightTree;
        newLeft.diffOrder = max(getDifforder(leftLeft), rightDiffOrder);
        newLeft.height = max(getHeight(leftLeft), rightHeight) + 1;
        newLeft.hasTerms = getHasTerms(leftLeft) || getHasTerms(rightTree);

        % Multiply together the new right factor of the left tree, and the
        % current right tree:
        newRight = treeVar();
        newRight.method = 'times';
        newRight.numArgs = 2;
        newRight.left = leftRight;
        newRight.right = rightTree;
        newRight.diffOrder = max(getDifforder(leftRight), rightDiffOrder);
        newRight.height = max(getHeight(leftRight), rightHeight) + 1;
        newRight.hasTerms = getHasTerms(leftRight) || getHasTerms(rightTree);
        
    elseif ( ~splitLeft && splitRight )
        % Had to split on the right, not the left.
        
        % Multiply together the current left tree, and the new left factor of
        % the right tree:
        newLeft = treeVar();
        newLeft.method = 'times';
        newLeft.numArgs = 2;
        newLeft.left = leftTree;
        newLeft.right = rightLeft;
        newLeft.diffOrder = max(leftDiffOrder, getDifforder(rightLeft));
        newLeft.height = max(leftHeight, getHeight(rightLeft)) + 1;
        newLeft.hasTerms = getHasTerms(leftTree) || getHasTerms(rightLeft);
        
        % Multiply together the current left tree, and the new right factor of
        % the right tree:
        newRight = treeVar();
        newRight.method = 'times';
        newRight.numArgs = 2;
        newRight.left = leftTree;
        newRight.right = rightRight;
        newRight.diffOrder = max(leftDiffOrder, getDifforder(rightRight));
        newRight.height = max(leftHeight, getHeight(rightRight)) + 1;
        newRight.hasTerms = getHasTerms(leftTree) || getHasTerms(rightRight);
        
    end
    
    % If the left or right trees still have multiple terms left in them, we need
    % to split them up recursively:
    if ( newLeft.hasTerms )
        newLeft = treeVar.expandTree(newLeft, maxOrder);
    end
    if ( newRight.hasTerms )
        newRight = treeVar.expandTree(newRight, maxOrder);
    end
    
    % Add together the new left and right trees, splitting complete!
    newTree = treeVar();
    newTree.method = 'plus';
    newTree.numArgs = 2;
    newTree.left = newLeft;
    newTree.right = newRight;
    newTree.diffOrder = max(newLeft.diffOrder, newRight.diffOrder);
    newTree.height = max(newLeft.height, newRight.height) + 1;
    newTree.hasTerms = newLeft.hasTerms || newRight.hasTerms;
    
end
end

function out = getDifforder(treeIn)
%GETDIFFORDER   The diffOrder of a TREEVAR.
%   We often encounter child nodes that may or may not be a syntax tree, as
%   opposed to scalars. We are usually interested in the heights and diffOrders
%   of the syntax trees, however, if the tree is actually a scalar, they won't
%   have the necessary fields. This method takes care of checking the class of
%   the argument, and return the correct information.
if ( isa(treeIn, 'treeVar') )
    out = treeIn.diffOrder;
else
    out = 0;
end

end

function out = getHeight(treeIn)
%GETHEIGHT   The height of a TREEVAR.
%   Same as GETDIFFORDER() method but for the height of a syntax tree.
if ( isa(treeIn, 'treeVar') )
    out = treeIn.height;
else
    out = 0;
end

end

function out = getHasTerms(treeIn)
%GETHASTERMS   Test if a TREEVAR contains a 'plus' or 'minus' node.
%   Same as GETDIFFORDER() method but for whether the syntax tree consists of
%   multiple terms.
if ( isa(treeIn, 'treeVar') )
    out = treeIn.hasTerms;
else
    out = 0;
end

end
