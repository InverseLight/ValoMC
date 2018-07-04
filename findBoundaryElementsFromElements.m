function output = findBoundaryElementsFromElements(H, BH)   
   sortedBH = transpose(sort(transpose(BH)));
   sortedH = transpose(sort(transpose(H)));
   if(size(H,2) == 4) 
       [test loc1] =ismember(sortedBH,sortedH(:,[1 2 3]),'rows'); 
       [test2 loc2] =ismember(sortedBH,sortedH(:,[1 2 4]),'rows'); 
       [test3 loc3] =ismember(sortedBH,sortedH(:,[1 3 4]),'rows'); 
       [test4 loc4] =ismember(sortedBH,sortedH(:,[2 3 4]),'rows'); 
       output=max([loc1 loc2 loc3 loc4]');
   else
       [test loc1] =ismember(sortedBH,sortedH(:,[1 2]),'rows'); 
       [test2 loc2] =ismember(sortedBH,sortedH(:,[1 3]),'rows'); 
       [test3 loc3] =ismember(sortedBH,sortedH(:,[2 3]),'rows'); 
       output=max([loc1 loc2 loc3]');
   end
end     
