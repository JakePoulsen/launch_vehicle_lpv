function [str,file] = symbolic_to_file(fileName,fun)
name = fileName;
f=sym(fun);
file=name+".m";
open = which(file);
fid = fopen(open,'w+');
% if fileID < 0
%    error('Failed to open myfile because: %s', message);
% end
fprintf(fid, '[\n');
arrayfun(@(ROWIDX) fprintf(fid, '%s,',f(ROWIDX,1:end-1)) + fprintf(fid, '%s;\n', f(ROWIDX,end)), (1:size(f,1)).');
fprintf(fid, '];\n');
fclose(fid);
str=fileread(file);
end