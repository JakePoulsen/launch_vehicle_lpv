str_start = 'function Bd_r=Bd_rigid_pitch(x_CG,x_PVP,S_ref,J_y,m,rho,v_x,idx,cd_coeff,aoa,x_CP)\n';
str_end = 'end';
str_Bd=replace(str_Bd,{'[','c_d_v(v_x, 0, 0, mach)'},{'Bd_r = [','polyval(cd_coeff{idx}, aoa)'});
open = fullfile(find_uncertain,file_Bd);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str_Bd);
fprintf(fid,str_end);
fclose(fid);