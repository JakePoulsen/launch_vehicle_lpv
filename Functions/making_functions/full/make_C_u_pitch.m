str_start = 'function C_u=C_rigid_pitch(x_CG,x_CP,x_INS,rho,v_x,S_ref,Psi_z1_INS,Psi_z2_INS,Psi_theta1_INS,Psi_theta2_INS)\n';
str_end = 'end';
str = replace(str_C_u, {'[', '(v_x^2)^(1/2)'}, {'C_u = [', 'v_x'});
open = fullfile(find_uncertain,file_C_u);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str);
fprintf(fid,str_end);
fclose(fid);