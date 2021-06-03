str_start = 'function Dn_r_u=Dn_rigid_pitch()\n';
str_end = 'end';
str_Dn=strrep(str_Dn,'[','Dn_r_u = [');
open = fullfile(find_uncertain,file_Dn);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str_Dn);
fprintf(fid,str_end);
fclose(fid);