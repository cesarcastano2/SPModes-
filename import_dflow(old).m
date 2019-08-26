function [Frame, Time, markers, forces] = import_dflow(dflow_file)

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
N = length(varnames);

% to read in numeric data
formatspec = repmat('%f ',1,N);
tempdata = textscan(file_id, sprintf('%s', formatspec), 'delimiter', ',');
tempdatamat = cell2mat(tempdata);
clear tempdata

fclose(file_id);

%% Frames & Time
Frame = tempdatamat(:,2); % Frame
startidx = find(Frame == 0);
stopidx = find(diff(Frame) == max(diff(Frame)));
if length(stopidx) ~= 1, stopidx = length(Frame); end
Time = tempdatamat(startidx:stopidx,1)-tempdatamat(startidx,1); % TimeStamp

Frame = tempdatamat(startidx:stopidx,2); % Frame

% lower limb markers
markers.idx = find(cellfun(@isempty,strfind(varnames, '.Pos')) == 0);
markers.fullnames = varnames(markers.idx);

ct = 1;
for m = 1:3:length(markers.idx)
    dotidx = find(markers.fullnames{m} == '.');   
    markers.labels{ct} = markers.fullnames{markers.idx(m)}(1:dotidx-1);
    ct = ct+1;
   
    eval(['markers.' markers.fullnames{markers.idx(m)}(1:dotidx-1) ' = (tempdatamat(startidx:stopidx,markers.idx(m:m+2)));']);   
end

% force data
forces.fp1idx = find(cellfun(@isempty,strfind(varnames, 'FP1')) == 0);
forces.fp2idx = find(cellfun(@isempty,strfind(varnames, 'FP2')) == 0);
forces.idx = [forces.fp1idx forces.fp2idx];
forces.fullnames = varnames(forces.idx);

ct = 1;
for f = 1:3:length(forces.idx)
    dotidx = find(forces.fullnames{f} == '.');
    forces.labels{ct} = [forces.fullnames{f}(1:dotidx-1) forces.fullnames{f}(dotidx+1:end-1)];
    ct = ct+1;
    
    eval(['forces.' [forces.fullnames{f}(1:dotidx-1) forces.fullnames{f}(dotidx+1:end-1)] ' = (tempdatamat(startidx:stopidx,forces.idx(f:f+2)));']);
end

