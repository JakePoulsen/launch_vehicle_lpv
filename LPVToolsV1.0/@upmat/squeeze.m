function E = squeeze(M)
% SQUEEZE  Removes singleton array dimensions
%
% SQUEEZE(M) removes all singleton array dimensions from M.
%
% See also: lpvelimiv.

% Check # of input/output arguments
nin = nargin;
nout = nargout;
error(nargchk(1, 1, nin, 'struct'))
error(nargoutchk(0, 1, nout, 'struct'))

% Find singleton array dimensions
niv = M.Domain.NumIV;
nad = numel(size(M))-2;
if nad<=2
    E =M;
    return;
end
Mdata = M.Data;         % [row col IVs AD]
szM = size(Mdata);
rmidx = find( szM(2+niv+1:end) == 1);

% Remove singleton array dimensions
keepidx = setdiff( 1:(2+niv+nad), 2+niv+rmidx );
Mdata = permute(M.Data,[keepidx, 2+niv+rmidx]);

% Pack up data
E = upmat(Mdata,M.Domain);


