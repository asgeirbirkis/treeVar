classdef  (InferiorClasses = {?chebfun}) treeVar
    %TREEVAR   A class for analysing syntax trees of ODEs
    %   The TREEVAR class is used for analysing the syntax trees of ODEs.
    %   Originally, it was implemented to be applied in CHEBFUN [1]. Its current
    %   use is to enable Chebfun to automatically convert (systems of) higher
    %   order ODEs to coupled first order systems. This is particularly useful
    %   for initial-value problems (IVPs), as that allows Chebfun to call one of
    %   the built-in MATLAB solvers for solving IVPs via time-stepping, rather
    %   than globally via spectral methods and Newton's method in function
    %   space.
    %
    %   T = TREEVAR(ID), where ID is a Boolean vector corresponding to the order
    %   of variables in the problem, returns the TREEVAR object T, which stores
    %   the ID. See example below for how the ID vector is specified.
    %
    %   T = TREEVAR() is the same as above, but with the default ID = 1. This is
    %   useful for quick testing purposes.
    %
    %   Example 1: Constructing TREEVAR object for the scalar IVP
    %       u'' + sin(u) = 0
    %   is done as follows:
    %       u = treeVar(1);
    %
    %   Example 2: Constructing TREEVAR objects for the coupled IVP
    %       u'' + v = 1, u + v' = x
    %   is done as follows:
    %       u = treeVar([1 0]);
    %       v = treeVar([0 1]);
    %
    % See also TREEVAR.SOLVEIVP.
    
    % References
    %  [1] Chebfun: http://www.chebfun.org
    
    % Copyright 2015 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TREEVAR class description:
    %
    % The TREEVAR class is used to convert higher order ODEs to coupled systems
    % of first order ODEs, which can then be solved using one of the built-in
    % MATLAB solvers, such as ODE113. This is done by evaluating (anonymous)
    % functions with TREEVAR arguments, which will construct a syntax tree of
    % the mathematical expression describing the operator. By then analysing the
    % syntax tree and restructuring it appropriately, conversion to a first-
    % order system is made possible.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        % The syntax tree of a TREEVAR variable, starting from an initial
        % variable. Each TREEVAR variable object contains the following fields:
        %
        %    METHOD:    The method leading to the construction of the variable.
        %    NUMARGS:   Number of arguments to the method that constructed the
        %        variable.
        %    DIFFORDER: The differential order of the TREEVAR, which
        %        represents how many times the base variable(s) have been
        %        differentiated when before we arrive at the current TREEVAR.
        %        Note that DIFFORDER is vector valued; for example, the
        %        sequence
        %            u = treeVar([1 0 0]);
        %            v = treeVar([0 1 0]);
        %            w = treeVar([0 0 1]);
        %            f = diff(u) + diff(w, 2);
        %        will lead to f.tree.DIFFORDER == [1 0 2].
        %    HEIGHT: The height of the syntax tree, i.e., the number of
        %        operations between the base variables(s) and the current
        %        variable.
        %    ID: A Boolean vector, whose ith element is equal to 1 if the
        %        TREEVAR variable was constructed from the ith base variable,
        %        0 otherwise. For example, the sequence
        %            u = treeVar([1 0 0]);
        %            v = treeVar([0 1 0]);
        %            w = treeVar([0 0 1]);
        %            f = u + 2*w;
        %        will lead to f.tree.ID == [1 0 1].
        %    HASTERMS: Indiciates whether a TREEVAR is constructed from a
        %        sequence of computations that include multiple terms. For
        %        example, the sequence
        %            u = treeVar();
        %            v = cos(x).*u;
        %            w = cos(u) + u;
        %        leads to v.tree.hasTerms = 0, w.tree.hasTerms = 1.
        diffOrder
        hasTerms = 0;
        height = 0;
        ID
        method = 'constr';
        numArgs = 0;
        left
        center
        right
    end
    
    properties( Hidden = true )
        x % For plotting
        y % For plotting
    end
    
    methods
        
        function obj = treeVar(IDvec)
            % The TREEVAR constructor. See documentation above for calling
            % sequences to the constructor.
            
            if ( nargin == 0 )
                % Default ID.
                IDvec = 1;
            end
            
            % Initialise a syntax tree for a base variable:
            obj.diffOrder = 0*IDvec;
            obj.ID = logical(IDvec);
        end
        
        function f = abs(f)
            f = univariate(f, 'abs');
        end
        
        function f = acos(f)
            f = univariate(f, 'acos');
        end
        
        function f = acosd(f)
            f = univariate(f, 'acosd');
        end
        
        function f = acot(f)
            f = univariate(f, 'acot');
        end
        
        function f = acoth(f)
            f = univariate(f, 'acoth');
        end
        
        function f = acsc(f)
            f = univariate(f, 'acsc');
        end
        
        function f = acscd(f)
            f = univariate(f, 'acscd');
        end
        
        function f = acsch(f)
            f = univariate(f, 'acsch');
        end
        
        function f = airy(f)
            f = univariate(f, 'airy');
        end
        
        function f = asec(f)
            f = univariate(f, 'asec');
        end
        
        function f = asecd(f)
            f = univariate(f, 'asecd');
        end
        
        function f = asech(f)
            f = univariate(f, 'asech');
        end
        
        function f = asin(f)
            f = univariate(f, 'asin');
        end
        
        function f = asind(f)
            f = univariate(f, 'asind');
        end
        
        function f = asinh(f)
            f = univariate(f, 'asinh');
        end
        
        function f = atan(f)
            f = univariate(f, 'atan');
        end
        
        function f = atand(f)
            f = univariate(f, 'atand');
        end
        
        function f = atanh(f)
            f = univariate(f, 'atanh');
        end
        
        function f = cos(f)
            f = univariate(f, 'cos');
        end
        
        function f = cosd(f)
            f = univariate(f, 'cosd');
        end
        
        function f = cosh(f)
            f = univariate(f, 'cosh');
        end
        
        function f = cot(f)
            f = univariate(f, 'cot');
        end
        
        function f = cotd(f)
            f = univariate(f, 'cotd');
        end
        
        function f = coth(f)
            f = univariate(f, 'coth');
        end
        
        function f = csc(f)
            f = univariate(f, 'csc');
        end
        
        function f = cscd(f)
            f = univariate(f, 'cscd');
        end
        
        function f = csch(f)
            f = univariate(f, 'csch');
        end
        
        function f = cumsum(f)
            %CUMSUM   Not supported.
            % We don't support integral equations with our first order
            % reformulation. However, we could accidentally end up here in case
            % of first order integral equation, where the conditions are
            % specified via N.LBC/RBC. Throw a meaningful error message in this
            % case.
            error('TREEVAR:cumsum:notSupported', ['First order ' ...
                'reformulation does currently not support integral equations.'])
        end
        
        function g = diff(f, k)
            %DIFF   Derivative of a TREEVAR.
            
            % By default, compute first derivative:
            if ( nargin < 2 )
                k = 1;
            end
            
            % Initialize
            g = f;
            
            g.method = 'diff';
            g.numArgs = 2;
            g.left = f;
            g.right = k;
            g.diffOrder = f.diffOrder + k*f.ID;
            g.height = f.height + 1;
            
        end
                
        function f = exp(f)
            f = univariate(f, 'exp');
        end
        
        function f = expm1(f)
            f = univariate(f, 'expm1');
        end
        
        function f = log(f)
            f = univariate(f, 'log');
        end
        
        function f = log10(f)
            f = univariate(f, 'log10');
        end
        
        function f = log2(f)
            f = univariate(f, 'log2');
        end
        
        function f = log1p(f)
            f = univariate(f, 'log1p');
        end
        
        function h = minus(f, g)
            %-   Subtraction of TREEVAR objects.
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                % SCALAR - TREEVAR
                h = bivariate(f, g, 'minus', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR - SCALAR
                h = bivariate(f, g, 'minus', 0);
            else
                % TREEVAR - TREEVAR
                h = bivariate(f, g, 'minus', 2);
            end
        end
        
        function h = mrdivide(f, g)
            %/   Matrix division of TREEVAR objects
            %
            % This method only supports (SCALAR/TREEVAR)/(SCALAR/TREEVAR)
            if ( isnumeric(f) || isnumeric(g) )
                h = rdivide(f, g);
            else
                error('Dimension mismatch');
            end
        end
        
        
        function h = mtimes(f, g)
            %*   Matrix multiplication of TREEVAR objects.

            % This method only supports SCALAR/TREEVAR*SCALAR/TREEVAR
            if ( isnumeric(f) || isnumeric(g) )
                h = times(f, g);
            else
                error('Dimension mismatch');
            end
        end
        
        function f = pow2(f)
            f = univariate(f, 'pow2');
        end
        
        function h = power(f, g)
            %.^   Power of a TREEVAR.
            if ( ~isa(f, 'treeVar') )
                % SCALAR.^TREEVAR
                h = bivariate(f, g, 'power', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR.^SCALAR
                h = bivariate(f, g, 'power', 0);
            else
                % TREEVAR.^TREEVAR
                h = bivariate(f, g, 'power', 2);
            end
        end
                
        function h = plus(f, g)
            if ( ~isa(f, 'treeVar') )
                % SCALAR+TREEVAR
                h = bivariate(f, g, 'plus', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR + SCALAR
                h = bivariate(f, g, 'plus', 0);
            else
                % TREEVAR + TREEVAR
                h = bivariate(f, g, 'plus', 2);
            end
        end
        
        function h = rdivide(f, g)
            %./   Division of TREEVAR objects.
            if ( ~isa(f, 'treeVar') )
                % SCALAR./TREEVAR
                h = bivariate(f, g, 'rdivide', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR./SCALAR
                h = bivariate(f, g, 'rdivide', 0);
            else
                % TREEVAR./TREEVAR
                h = bivariate(f, g, 'rdivide', 2);
            end
        end
        
        function f = sec(f)
            f = univariate(f, 'sec');
        end
        
        function f = secd(f)
            f = univariate(f, 'secd');
        end
        
        function f = sech(f)
            f = univariate(f, 'sech');
        end

        function f = sin(f)
            f = univariate(f, 'sin');
        end
        
        function f = sind(f)
            f = univariate(f, 'sind');
        end
        
        function f = sinh(f)
            f = univariate(f, 'sinh');
        end

        function f = sqrt(f)
            f = univariate(f, 'sqrt');
        end
        
        function f = tan(f)
            f = univariate(f, 'tan');
        end
        
        function f = tand(f)
            f = univariate(f, 'tand');
        end
        
        function f = tanh(f)
            f = univariate(f, 'tanh');
        end
        
        function h = times(f, g)
            %.*   Multiplication of treeVar objects.
            
            % Initialize an empty TREEVAR
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                % SCALAR.^*TREEVAR
                h = bivariate(f, g, 'times', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR.*SCALAR
                h = bivariate(f, g, 'times', 0);
            else
                % TREEVAR.*TREEVAR
                h = bivariate(f, g, 'times', 2);
            end
        end
        
        function f = uminus(f)
            f = univariate(f, 'uminus');
        end
        
        function f = uplus(f)
            f = univariate(f, 'uplus');
        end
    end

    
    methods ( Static = true )
        
        % Setup initial conditions for solving an IVP
        initVals = setupICs(icFun, totalDiffOrders)
        
        % Solve an IVP
        [t, y] = solveIVP(odeFun, icFun, rhs, dom)
        
        % Returns how the results of evaluating BCs should be sorted
        idx = sortConditions(funIn, domain, maxDiffOrders)
        
        % Convert higher order anonymous functions to first order systems
        [funOut, indexStart, problemDom, coeffs, totalDiffOrders] = ...
            toFirstOrder(funIn, rhs, domain)
        
    end
    
    methods ( Static = true, Access = private )
                
        % Convert expressions like 5*(diff(u) + u) to 5*diff(u) + 5*u
        newTree = expandTree(tree, maxOrder)
                
        % Split syntax trees into derivative part and non-derivative part
        [newTree, derTree] = splitTree(tree, maxOrder)
        
        % Convert the infix form of an expression to an anonymous function
        anonFun = toAnon(infix, varArray)
        
        % Convert infix expressions to anonymous function suited for ODE solvers
        funOut = toRHS(infix, varArray, coeff, indexStart, totalDiffOrders);
        
        % Convert a syntax tree to infix form
        [infix, varArray] = ...
            tree2infix(tree, diffOrders, varCounter, varArray, isCoeffFun)
                   
    end

end
