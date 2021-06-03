str_start = 'function C_r_u=C_rigid_pitch(x_CG,x_INS,rho,v_x)\n';
str_end = 'end';
str = replace(str_C, {'[', '(v_x^2)^(1/2)'}, {'C_r_u = [', 'v_x'});
open = fullfile(find_uncertain,file_C);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str);
fprintf(fid,str_end);
fclose(fid);