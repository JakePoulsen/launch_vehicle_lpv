
blockName='Clock1'    
sys=load_system(gcs+".slx");
b = getfullname(Simulink.findBlocks(sys))
%contains(b,'Clock')
%find(ismember(b, 'dynamics'))
index = find(~cellfun(@isempty, regexp(b,blockName+"$")))
paths=b(index')

for p=1:size(paths,1)
   p 
end