function [P]  = findLineSegmentsFromElements(H, BH)
% Finds the element to which a line segment belongs to
%
% function [P]  = findLineSegmentsFromElements(H, BH)
%
% INPUT
%
%  H:		topology matrix
%  BH:	        boundary topology matrix
%
% OUTPUT
%
%  P:           vector of rows in H to which each row in BH belongs to. [size(BH,1)]
%
%
% TODO
%
%  vectorize, optimize
%
   P = ones(length(BH),1);

   for i=1:size(BH,1)
      for j=1:size(H,1)
           if(length(find(H(j,:) == BH(i,1))) && length(find(H(j,:) ...
                                                             == ...
                                                             BH(i, ...
                                                                2))))
               P(i) = j;
           end
      end
   end
end

