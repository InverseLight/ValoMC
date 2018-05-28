function [HN] = createHN(H);
%CREATEHN Creates topology neighbourhood HN for topology H
%
%
% DESCRIPTION:
%       This can function can be used to create a neighborhood
%       matrix for a given element topology
%
% USAGE:
%       HN = createHN(H)
%
% INPUTS:
%       H           - Element topology
%
% OUTPUTS:
%       HN          - Neighborhood topology



if(size(H,2) > 3)
   error('3d meshes are not currently supported')
end

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

end;
