function dimension = getInputDimension(filename)
   fileID = fopen(filename,'r');
   tline = fgets(fileID);
   tline = fgets(fileID);
   tline = fgets(fileID);
   tline = fgets(fileID);
   tline = fgets(fileID);
   C = find(~cellfun(@isempty,strsplit(tline)));
   fclose(fileID);
   dimension=0;
   if(size(C,2)==3) 
       dimension=2;
   else if(size(C,2)==4) 
       dimension=3;
   end
end
