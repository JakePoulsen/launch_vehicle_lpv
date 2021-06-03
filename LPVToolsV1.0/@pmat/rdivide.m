function out = rdivide(A,B)
% RDIVIDE  Right array divide for PMAT objects
%
% RDIVIDE(A,B) is the result of A./B at each point in the combined
% domains of A and B.
%
% See also: rdivide, ldivide, mldivide, mrdivide.

% Check # of input arguments
error(nargchk(2, 2, nargin, 'struct'))
out = binop(A,B,'rdivide');
