function handle = CodeImporter(blockName,fileName)
      %  open_system(gcs+".slx");
        sf = sfroot();
    sys=load_system(gcs+".slx");
    b = getfullname(Simulink.findBlocks(sys));
    %contains(b,'Clock')
    %find(ismember(b, 'dynamics'))
    index = find(~cellfun(@isempty, regexp(b,"/"+blockName+"$")));
    blockPaths=b(index');
     
    for p=1:size(blockPaths,1)
%         
          blockPaths{p}
        block = sf.find('Path',blockPaths{p},'-isa','Stateflow.EMChart');
        block.Script = fileread(fileName);
    end
    
%      blockPath = gcs+"/"+blockName;
%     % handle= getSimulinkBlockHandle(blockPath);
%      block = sf.find('Path',blockPath,'-isa','Stateflow.EMChart');
%      block.Script = fileread(fileName);
end
