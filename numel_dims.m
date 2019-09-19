function N = numel_dims(data,dims) 
% Compute number of elements in selected dimensions
%   Created by Yibo Zhao @ UIUC, 07/23/2018


    cmd = sprintf('N = numel(data(');
    for n = 1:ndims(data)
        if n ~= ndims(data)
            if any(n == dims)
                cmd = sprintf('%s:,',cmd);
            else
                cmd = sprintf('%s1,',cmd);
            end
        else
            if any(n == dims)
                cmd = sprintf('%s:',cmd);
            else
                cmd = sprintf('%s1',cmd);
            end
        end
    end
    cmd = sprintf('%s));',cmd);
    eval(cmd);
    
end

%#ok<*STOUT>

