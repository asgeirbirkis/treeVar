function treeOut = bivariate(leftTree, rightTree, method, type)
%BIVARIATE   Construct a syntax tree for bivariate functions.
%
% The BIVARIATE() method is called by most bivariate methods of the TREEVAR
% class, such as minus, plus and rtimes. It constructs the corresponding syntax
% tree, where the appropriate method is at the top, and the syntax trees that
% the method operates on become the left and right nodes.
%
% Calling sequence:
%   TREEOUT = BIVARIATE(LEFTTREE, RIGHTTREE, METHOD, TYPE)
%
% The inputs to this method are:
%   LEFTTREE:   A TREEVAR for representing the syntax tree of the mathematical
%               expressions that is the left argument that METHOD operates on.
%   RIGHTTREE:  A TREEVAR for representing the syntax tree of the mathematical
%               expressions that is the right argument that METHOD operates on.
%   METHOD:     The name of the method that invokes the call to BIVARIATE().
%   TYPE:       An integer that indicated which of the arguments to METHOD were
%               TREEVAR objects (as opposed to scalars). The possible values
%               are:
%                   0: Only the left argument was a TREEVAR,
%                   1: Only the right argument was a TREEVAR,
%                   2: Both left and right arguments were TREEVAR objects.   
%
% The output of this method is:
%   TREEOUT: A TREEVAR that represents the syntax tree of the mathematical
%            expression, obtained once METHOD has operated on LEFTTREE and
%            RIGHTTREE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Are we dealing with a plus or minus (which affects the HASTERMS field of the
% output).
isPM = any( strcmp(method, {'plus', 'minus'}) );

if ( type == 2 )
    % Both LEFTTREE and RIGHTTREE were TREEVAR objects.
    treeOut = treeVar();
    treeOut.ID = leftTree.ID | rightTree.ID;
    treeOut.diffOrder = max(leftTree.diffOrder, rightTree.diffOrder);
    treeOut.height = max(leftTree.height, rightTree.height) + 1;
    
    % We only care about terms for dependent variables
    hasDepVars = any(leftTree.ID) || any(rightTree.ID);
    treeOut.hasTerms = hasDepVars && (isPM || leftTree.hasTerms ...
        || rightTree.hasTerms);    
    
elseif ( type == 1 )
    % Only RIGHTTREE was a TREEVAR.
    treeOut = rightTree;
    treeOut.height = treeOut.height + 1;
    hasDepVars = any(rightTree.ID);
    treeOut.hasTerms = hasDepVars && ( isPM || rightTree.hasTerms );
else
    % Only LEFTTREE was a TREEVAR.
    treeOut = leftTree;
    treeOut.height = treeOut.height + 1;
    hasDepVars = any(leftTree.ID);
    treeOut.hasTerms = hasDepVars && ( isPM || leftTree.hasTerms );
end

treeOut.numArgs = 2;
treeOut.method = method;
treeOut.left = leftTree;
treeOut.right = rightTree;

end
