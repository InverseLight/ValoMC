function [value] = readNetgenValue(entryname, fid, tline)
% Auxiliary function for importNetGenmesh
%
%  This function is used to extract scalars from the vol files
%  (like the dimensionality).
%
    value = [];
    if(strcmp(tline, entryname))
        value = str2num(fgetl(fid));
    end
end
