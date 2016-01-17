function treeOut = univariate(treeIn, method)
%UNIVARIATE   Construct a syntax tree for univariate functions.
%   The UNIVARIATE method is called by most univariate methods of the TREEVAR
%   class, such as COS, EXP and SIN. It constructs the corresponding syntax
%   tree, where the appropriate method is at the top and the syntax tree up to
%   that point becomes the center node.
%
%   Calling sequence:
%      TREEOUT = UNIVARIATE(TREEIN, METHOD)
%
%   The inputs to this method are:
%      TREEIN: A TREEVAR for representing the syntax tree of a mathematical
%              expressions that the METHOD operates on.
%      METHOD: The name of the method that invokes the call to UNIVARIATE.
%
%   The output of this method is:
%      TREEOUT: A TREEVAR that represents the syntax tree of the mathematical 
%                expression, obtained once METHOD has operated on TREEIN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Construct a new TREEVAR to be returned.
treeOut = treeIn;
treeOut.height = treeIn.height + 1;
treeOut.method = method;
treeOut.numArgs = 1;
treeOut.center = treeIn;
end
