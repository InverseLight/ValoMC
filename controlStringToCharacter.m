function c = controlStringToCharacter(string, default)
%CONTROLSTRINGTOCHARACTER Convert a control string to a single character for the C++ code 
%
%
% DESCRIPTION:
%       (ValoMC internal use only) 
%       Converts strings that are used in the ValoMC Matlab interface like 'cosinic'
%       into a single byte character. For example 'gaussian' -> 'g'
%   
% USAGE:
%       output = controlStringToCharacter('gaussian', 'n')
%
% INPUTS:
%       string     - Control string to convert
%       default    - If no character is found for conversion, return this value
%
% OUTPUTS:
%       c          - Converted character

  if(strcmp(string,'none')) c = 'a';
  elseif(strcmp(string,'direct')) c = 'l';
  elseif(strcmp(string,'cosinic')) c = 'c';
  elseif(strcmp(string,'gaussian')) c = 'g';
  elseif(strcmp(string,'isotropic')) c = 'i';
  elseif(strcmp(string,'absolute')) c = 'a';
  elseif(strcmp(string,'relative')) c = 'r';
  elseif(strcmp(string,'pencil')) c = 'p';
  else c = default;
  end

end

