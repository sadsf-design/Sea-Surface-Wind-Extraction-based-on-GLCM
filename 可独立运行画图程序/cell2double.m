function Out = cell2double(in)
    [r, c] = size(in);
    Out = zeros(r,c);
    for i = 1 : r
        for j = 1 : c
            if isa(in{i, j}, 'double')
               Out(i, j) = in{i, j};
            else
               Out(i, j) = str2double(in{i, j}); 
            end
        end
    end
end