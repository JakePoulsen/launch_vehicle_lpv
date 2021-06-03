str_start = 'function D_r_u=D_rigid_pitch()\n';
str_end = 'end';
str_D=strrep(str_D,'[','D_r_u = [');
open = fullfile(find_uncertain,file_D);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str_D);
fprintf(fid,str_end);
fclose(fid);