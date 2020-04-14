% Read raw equation and add '...' to be read by Matlab.

clc; clear;

% Read txt file
inpCell = regexp(fileread('eqRaw.txt'), '\n', 'split');

% Edit string
str = strjoin(inpCell);
strtrim(str);    % Remove leading and trailing white space from string
str(str==' ') = [];
str(str==sprintf('\t')) = [];    % Remove tab
str(str==sprintf('\r')) = [];    % Remove format
str(str==':') = [];
str = strrep(str, 'Abs', 'abs');
str = strrep(str, 'Exp', 'exp');
str = strrep(str, '^', '.^');
str = strrep(str, '*', '.*');

% Wrap string
eqProcessed = WRAP_string(str);

% Write txt file
fid     = fopen('eqProcessed.txt', 'w');
for i = 1:size(eqProcessed, 1)
    str = eqProcessed(i, :);
    fprintf(fid, '%s\r\n', str);
end
fclose(fid);