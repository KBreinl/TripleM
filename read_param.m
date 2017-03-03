function out = read_param(paramName)
% READ_PARAM
%
% out = read_param(paramName)
%
% Given a parameter name in the 'param.yml' file, return its value.
%

fileID = fopen('param.yml');
C = textscan(fileID, '%[^: ]%*[: ]%s', 'CommentStyle', '#');
fclose(fileID);

tmp = C{2}(strcmp(C{1}, paramName));

out = cell2mat(tmp);

% if the parameter is 'dist' or 'cl_period', return a string instead
if ~any(strncmp(paramName, {'dist', 'cl_period', 'corr_rand', 'para_rainf'}, 4))
    out = str2double(out);
end

end

