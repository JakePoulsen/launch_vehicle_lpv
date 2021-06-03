str_start = 'function Bd_r_u=Bd_rigid_pitch(x_CG,x_PVP,S_ref,J_y,m,mach,rho,v_x,x_CP,cd_coeff,idx,aoa)\n';
str_end = 'end';
str_Bd=replace(str_Bd_u,{'[','c_d_v(v_x, 0, 0, mach)'},{'Bd_r_u = [','polyval(cd_coeff{idx}, aoa)'});
open = fullfile(find_uncertain,file_Bd_u);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str_Bd);
fprintf(fid,str_end);
fclose(fid);