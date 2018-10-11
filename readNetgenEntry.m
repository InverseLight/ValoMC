function [arr_out] = readNetgenEntry(entryname, fid, tline)
    arr_out = [];
    value = [];
    if(strcmp(tline, entryname))
        number_of_entries = str2num(fgetl(fid));
        for i=1:number_of_entries
            curline = fgetl(fid);
            linevector = str2num(curline);
            if(isempty(linevector))
                [value] = strsplit(curline);
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
