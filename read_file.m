function [dates, data, nos] = read_file(filename)
% READ_FILE
%
% [dates, data, nos] = read_file(filename)
%
% Read a csv file containing a given number of
% precipitation data series and their associated dates
%
% FILE FORMAT:
%
% Date, S1, S2, ...
% yyyy-MM-dd, P11, P21, ...
% yyyy-MM-dd, P12, P22, ...
% yyyy-MM-dd, P13, P23, ...
% ...
%
% Out:
% dates - array of datetime objects
% data - matrix of precipitation data (one site per column)
% nos - number of sites
%

% Get the first line of the file and count the stations
fid = fopen(filename);
line = fgetl(fid);
fclose(fid);

tmp = textscan(line, '%s', 'Delimiter', ',');

nos = length(tmp{1})-1;  % Count the stations

% Scan the dates and convert them to datetime objects
fid = fopen(filename);
dates = textscan(fid, '%s%*[^\n]', 'Delimiter', ',', ...
                 'HeaderLines', 1);
fclose(fid);

dates = datetime(dates{1}, 'InputFormat', 'yyyy-MM-dd');

% Scan the data and convert it to a matrix
fid = fopen(filename);
data = textscan(fid, ['%*s ', repmat('%f ', 1, nos)], ...
                'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);

data = cell2mat(data);

end

