function [HN] = createHN(H);
%CREATEHN Creates topology neighbourhood HN for topology H
%
% USAGE:
%
%       HN = createHN(H)
%
% DESCRIPTION:
%
%       This can function can be used to create a neighborhood
%       matrix for a given element topology
%
% INPUT:
%
%       H           - Element topology
%
% OPTIONAL INPUT:
%
%       optarg1     - Description of the first optional argument [J/m3]
%
% OUTPUT:
%
%       HN          - Neighborhood topology
%
% This function is provided with ValoMC

if(size(H,2) > 3)
   error('Three dimensional meshes not yet supported');
else

HN = zeros(size(H));

for(el = 1 : size(H, 1));

  ii = setdiff( find( sum( reshape( ismember(H, H(el, [1 2])), size(H) ...
                                    ), 2 ) == 2 ), el);
  if(~isempty(ii));
    HN(el, 1) = ii;
  end;
  
  ii = setdiff( find( sum( reshape( ismember(H, H(el, [2 3])), size(H) ...
                                    ), 2 ) == 2 ), el);
  if(~isempty(ii));
    HN(el, 2) = ii;
  end;

  ii = setdiff( find( sum( reshape( ismember(H, H(el, [3 1])), size(H) ...
                                    ), 2 ) == 2 ), el);
  if(~isempty(ii));
    HN(el, 3) = ii;
  end;
end

end;
