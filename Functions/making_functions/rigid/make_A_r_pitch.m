str_start = 'function A_r_u=A_rigid_pitch(x_CG,x_CP,J_y,S_ref,m,mach,rho,v_x,g,idx,cd_coeff,cl_coeff,aoa)\n';
str_end = 'end';
str = replace(str_A, {'[','conj(x_CP)', 'conj(x_CG)','subs(diff(c_l_v(v_x, 0, v_z, mach), v_z), v_z, 0)', 'c_d_v(v_x, 0, 0, mach)','c_d'}, {'A_r_u = [','x_CP', 'x_CG','polyval(diff(cl_coeff{idx}), aoa)','polyval(cd_coeff{idx}, aoa)','polyval(cd_coeff{idx}, aoa)'});
open = fullfile(find_uncertain,file_A);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str);
fprintf(fid,str_end);
fclose(fid);