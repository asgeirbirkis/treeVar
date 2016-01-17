fprintf('test_bivariate      ...  ')
pass = test_bivariate;
if (all(pass))
    fprintf('passed!\n')
else
    fprintf('failed!\n')
end

fprintf('test_LorenzIVP      ...  ')
pass = test_LorenzIVP;
if (all(pass))
    fprintf('passed!\n')
else
    fprintf('failed!\n')
end

fprintf('test_plotTree       ...  ')
pass = test_plotTree;
if (all(pass))
    fprintf('passed!\n')
else
    fprintf('failed!\n')
end

fprintf('test_printTree      ...  ')
pass = test_printTree;
if (all(pass))
    fprintf('passed!\n')
else
    fprintf('failed!\n')
end

fprintf('test_sortConditions ...  ')
pass = test_sortConditions;
if (all(pass))
    fprintf('passed!\n')
else
    fprintf('failed!\n')
end

fprintf('test_toFirstOrder   ...  ')
pass = test_toFirstOrder;
if (all(pass))
    fprintf('passed!\n')
else
    fprintf('failed!\n')
end

fprintf('test_univariate     ...  ')
pass = test_univariate;
if (all(pass))
    fprintf('passed!\n')
else
    fprintf('failed!\n')
end