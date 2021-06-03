str_start = 'function Dd_r_u=Dd_rigid_pitch(rho,v_x)\n';
str_end = 'end';
str=replace(str_Dd,{'[','(v_x^2)^(1/2)'},{'Dd_r_u = [','v_x'});
open = fullfile(find_uncertain,file_Dd);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str);
fprintf(fid,str_end);
fclose(fid);