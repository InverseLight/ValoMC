function array_out  = extendCellArray(array_in, desired_size)
%COMPLETECELLARRAY Changes the size a cell array 
%
%
% DESCRIPTION:
%       A helper function to quicky take take content of one
%       array and make another, bigger cell array with it.
%
%
% USAGE:
%
% INPUTS:
%       array_in      Input array
%       desired_size  Desired size for the array
%
%
% OUTPUTS:
%       array_out      Extended array

  if(length(array_in) > desired_size) 
     warning('extendCellArray is intended only for extending arrays.');
  end
  array_out = cell(desired_size, 1);
  array_out(1:length(array_in)) = array_in(1:length(array_in));
end

