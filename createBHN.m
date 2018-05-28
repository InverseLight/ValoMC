function [BHN] = createBHN(BH);
%CREATEHN Creates topology neighbourhood BHN for boundary topology BH
%
%
% DESCRIPTION:
%       This can function can be used to create a neighborhood
%       matrix for a given element topology
%
% USAGE:
%       BHN = createBHN(BH)
%
% INPUTS:
%       BH           - Element topology
%
% OUTPUTS:
%       BHN          - Neighborhood topology



if(size(BH,2) > 2)
   error('3d meshes are not currently supported')
end

HN = zeros(size(BH));

for(el = 1 : size(BH, 1));

  ii = setdiff( find( sum( reshape( ismember(BH, BH(el, 1)), size(BH) ...
                                    ), 2 ) == 1 ), el);
  if(~isempty(ii));
    BHN(el, 1) = ii;
  end;

  ii = setdiff( find( sum( reshape( ismember(BH, BH(el, 2)), size(BH) ...
                                    ), 2 ) == 1 ), el);
  if(~isempty(ii));
    BHN(el, 2) = ii;
  end;

  
end;
