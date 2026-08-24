function out = allClose(x, y, options)
% This is the equivalent implementation of numpy.testing.assert_allclose
arguments
    x {isscalar, mustBeNumeric}
    y {isscalar, mustBeNumeric}
    options.atol {isscalar, mustBeNumeric} = 0.0
    options.rtol {isscalar, mustBeNumeric} = 1e-05
end
out = all(abs(x - y) <= options.atol + options.rtol * abs(y));
end