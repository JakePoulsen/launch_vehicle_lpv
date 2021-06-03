function out = vertcat(varargin)
% VERTCAT  Vertical concatenation of UPSS objects.
%
% SYS = VERTCAT(SYS1,SYS2,...) performs a concatenation operation of
% SYS = [SYS1; SYS2; ...]  at each point in the combined domains
% of SYS1, SYS2, ...
%
% See also: vertcat, horzcat.

% Check # of input arguments
error(nargchk(1, inf, nargin, 'struct'))

if nargin==1
    out = varargin{1};
else
    varargin{1} = upss(varargin{1});
    out = binop(varargin{1},varargin{2},'vertcat');
    if nargin>2
        out = vertcat(out,varargin{3:end});
    end    
end

