str=replace(str_Dd,{'[','(v_x^2)^(1/2)'},{'Dd_r_lpv = [','v_x'});
open = fullfile(find_lpv_rigid,file_Dd);
fid = fopen(open,'w+');
fprintf(fid,str_start);
fprintf(fid,str);
fprintf(fid,str_end);
fclose(fid);