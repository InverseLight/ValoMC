function array_out = duplicateArray(array_in, desired_size, defaultvalue)
%
% Duplicate array_in unit it has a desired size
%
% function array_out = duplicateArray(array_in, size)
%
% INPUT
%
%  array_in:		incomplete array
%  desired_size:	desired size for the array
%
% OUTPUT
%
%  array_out:   	completed array [size]
%

   % convert to column vector
   if(size(array_in,2) > size(array_in,1)) 
       array_in = transpose(array_in);
   end

   if(size(array_in,1) > desired_size)
      tmp = array_in(1:desired_size);
   else
      tmp = repmat(array_in,ceil(desired_size/size(array_in,1)),1);
   end
   
   array_out = tmp(1:desired_size,:);
   
end

