function [X,Y,Z,VX,VY,VZ] = txt2csv(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables

X = {};
Y = {};
Z = {};
VX = {};
VY = {};
VZ = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the file, delete unnecessary text, and read numeric data into
% cell array

fid = fopen(filename);  % Open the file

% Skip over the first 14 lines, since it doesn't contain any data
for i = 1:49
    fgetl(fid);
end

k = 1;  % Index the line

while ~feof(fid)
    tline = fgetl(fid); % Get the index of the current line
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input data into initialized data and vectors
    
    X{k}    = tline(52:73);
    Y{k}    = tline(76:97);
    Z{k}    = tline(100:121);
    VX{k}   = tline(124:145);
    VY{k}   = tline(148:169);
    VZ{k}   = tline(172:193);
    
    k = k + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finish reading document
    
    if strcmp(tline,'$$EOE')
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store cell string data in numeric vectors
X = str2double(X);
Y = str2double(Y);
Z = str2double(Z);
VX = str2double(VX);
VY = str2double(VY);
VZ = str2double(VZ);

end



