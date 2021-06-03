str_start = 'function B_r_u=B_rigid_pitch(x_CG,x_PVP,T_TVC,J_y,m)\n';
str_end = 'end';
str=strrep(str_B,'[','B_r_u = [');
open = fullfile(find_uncertain,file_B);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str);
fprintf(fid,str_end);
fclose(fid);