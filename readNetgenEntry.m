function [arr_out] = readNetgenEntry(entryname, fid, tline)
% Auxiliary function for importNetGenmesh
%
%  This function is used to extract entries from vol files
%  into matlab arrays
%
    arr_out = [];
    value = [];
    if(strcmp(tline, entryname))
        number_of_entries = str2num(fgetl(fid));
        for i=1:number_of_entries
            curline = fgetl(fid);
            linevector = str2num(curline);
            if(isempty(linevector))
                % reading a line as numbers was not succesfull
                % this is likely a number-string entry
                [value] = strsplit(curline);
                % return value is a cell array
                arr_out{uint32(str2num(value{1}))} = value{2};
            else
                arr_out = [arr_out; linevector];
            end
        end
        if(not(ischar(tline)))
            error('Could not read file');
        end
    end
end
