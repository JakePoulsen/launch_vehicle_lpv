str_start = 'function Bn_r=Bn_rigid_pitch(J_ny,m_n,x_PVP,x_n,x_CG,J_y,m)\n';
str_end = 'end';
str_Bn=replace(str_Bn,{'['},{'Bn_r = ['});
open = fullfile(find_uncertain,file_Bn);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str_Bn);
fprintf(fid,str_end);
fclose(fid);