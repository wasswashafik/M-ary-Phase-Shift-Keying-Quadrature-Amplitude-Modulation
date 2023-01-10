function y = qfunc(x)
%QFUNC  Q function.
%   Y = QFUNC(X) returns 1 minus the cumulative distribution function of the 
%   standardized normal random variable for each element of X.  X must be a real
%   array. The Q function is defined as:
%
%     Q(x) = 1/sqrt(2*pi) * integral from x to inf of exp(-t^2/2) dt
%
%   It is related to the complementary error function (erfc) according to
%
%     Q(x) = 0.5 * erfc(x/sqrt(2))
%
%   See also QFUNCINV, ERF, ERFC, ERFCX, ERFINV, ERFCINV.

%   Copyright 1996-2011 The MathWorks, Inc.

if (~isreal(x) || ischar(x))
  error(message('comm:qfunc:InvalidArg')); 
end
y = 0.5 * erfc(x/sqrt(2));
return;