function vout = import_dflow_treadmill(dflow_file)

%% Imports data from an D-FLOW txt file into matlab
%
% INPUT VARIABLES:
% dflow_file = full path for the DFlow file to be analyzed
%
% OUTPUT VARIABLES:
% sysinfo = contains information from the first row of the file with
% information such as the sampling frequency, file name, etc.
% Frame = vector of the frame numbers for the file
% Time = vector of the corresponding times [sec]
%

file_id = fopen(dflow_file);

% Get column names
dfheader = fgetl(file_id);
varnames = regexp(dfheader, '\t', 'split');
varnames(2)={'LeftBeltSpeed'};
varnames(3)={'RightBeltSpeed'};
N = length(varnames);

% to read in numeric data
formatspec = repmat('%f ',1,N);
tempdata = textscan(file_id, sprintf('%s', formatspec), 'delimiter', ',');
tempdatamat = cell2mat(tempdata);
clear tempdata

fclose(file_id);

%% 
for v = 1:length(varnames)
    spaceidx = strfind(varnames{v}, ' ');
    varnames{v}(spaceidx) = [];
    dashidx = strfind(varnames{v}, '-');
    varnames{v}(dashidx) = [];
    eval([ 'vout.' varnames{v} ' = tempdatamat(:,v);'])
end


