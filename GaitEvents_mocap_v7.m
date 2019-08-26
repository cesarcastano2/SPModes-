function[RHS,LTO,LHS,RTO,GE,GEInMiddleDeleted] = GaitEvents_mocap_v7(Time, llmarkers, markers4GE, HSrefinePre, HSrefinePost, TOminpeakdistance, TOminpeakheight,BW)   
% Inputs are:
%  Time: "TimeDF" from Dflow file
%  llmarkers: "markers" from Dflow file amd add "markers.labels(1,:) = {'LASIS', 'RASIS', 'LPSIS', 'RPSIS', 'LKNE', 'LTHI', 'LANK', 'LTIB', 'LTOE', 'LHEE', 'RKNE', 'RTHI', 'RANK', 'RTIB', 'RTOE', 'RHEE'}; " if use 16 markers
%  markers4GE: "{'RHEE' 'LHEE' 'RANK' 'LANK' 'RTOE' 'LTOE'}"
%  HSrefinePre: search window length before heel strike, must be positive, suggest use default as 10
%  HSrefinePost: search window length after heel strike, must be positive, suggest use default as 5
%  TOminpeakdistance: peak distance threshold for searching toe off, must be positive, suggest use default as 40
%  TOminpeakheight: peak height threshold for searching toe off, must be positive, suggest use default as 4
%  BW: use 0 for normal waling and 1 for backwards walking
% Output is a matrix GE if the Giat Event can be aligned precisely as:
%   Column 1 contains all time points of right heel strike
%   Column 2 contains all time points of left toe off
%   Column 3 contains all time points of left heel strike
%   Column 4 contians all time points of right heel strike
%   i.e. each row is a complete gait cycle startig with right heel strike.
%   if the Giat Event cannot be aligned precisely then output RHS,LTO,LHS,RTO respectively
%   Plus, GEInMiddleDeleted shows the Gait Event index which was deleted in the middle 

%%
if nargin >= 9
    error('Too many inputs');
elseif nargin == 8
    if BW == 1 % analyze backwards walking data in reverse, to function like forwards walking
        for m = 1:length(llmarkers.labels)
            eval(['llmarkers.' llmarkers.labels{m} ' = flipud(llmarkers.' llmarkers.labels{m} ');']);     
        end
        Time = flipud(Time);
    end
elseif nargin <= 7
    error('Not enough inputs');
end

%% Filter data, zero-lag fourth order low pass Butterworth filter at 12 Hz
ts = 1/240;  
fc = 12;     

[B_coef_filter, A_coef_filter] = butter(4,2.0*ts*fc );       
for i=1:4
    eval([markers4GE{i} 'yF =filtfilt(B_coef_filter,A_coef_filter,llmarkers.' markers4GE{i} '(:,2));']);    % y and z of 'RHEE' 'LHEE' 'RANK' 'LANK' for HS
    eval([markers4GE{i} 'zF =filtfilt(B_coef_filter,A_coef_filter,llmarkers.' markers4GE{i} '(:,3));']);
end
for i=5:6
    eval([markers4GE{i} '=filtfilt(B_coef_filter,A_coef_filter,llmarkers.' markers4GE{i} '(:,2));']);       % y of 'RTOE' 'LTOE' for TO
end

%% heel strike and toe off parameter
h = 1/240;               
Time = Time(4:end-3);    
m = size(RTOE,1);

%% Heel Strike    
for i = 4:m-3        % m = size(RTOE,1);
    VheelR(i-3) = (-RHEEyF(i+2) + 8*RHEEyF(i+1) - 8*RHEEyF(i-1) + RHEEyF(i-2))/(12*h);   
    VheelL(i-3) = (-LHEEyF(i+2) + 8*LHEEyF(i+1) - 8*LHEEyF(i-1) + LHEEyF(i-2))/(12*h);
    ACCheelR(i-3) = (-RHEEyF(i+2) + 16*RHEEyF(i+1) - 30*RHEEyF(i) + 16*RHEEyF(i-1) - RHEEyF(i-2) )/(12*h^2);                    
    ACCheelL(i-3) = (-LHEEyF(i+2) + 16*LHEEyF(i+1) - 30*LHEEyF(i) + 16*LHEEyF(i-1) - LHEEyF(i-2) )/(12*h^2);
    ACCJerkheelR(i-3) = (-RHEEyF(i+3) + 8*RHEEyF(i+2) - 13*RHEEyF(i+1) + 13*RHEEyF(i-1) - 8*RHEEyF(i-2) + RHEEyF(i-3))/(8*h^3); 
    ACCJerkheelL(i-3) = (-LHEEyF(i+3) + 8*LHEEyF(i+2) - 13*LHEEyF(i+1) + 13*LHEEyF(i-1) - 8*LHEEyF(i-2) + LHEEyF(i-3))/(8*h^3);
end

% left
[pks,leftHSpeak] = findpeaks(LHEEyF(1:end),'minpeakdistance',150);  
j=0;
for i = 1:length(leftHSpeak)-1   
    in = leftHSpeak(i);
    fi = leftHSpeak(i+1);        
    [A,B] = min(LHEEyF(in:fi));  
    j=j+1;
    NearleftHSi(j) = B+in-1;     
    
    CutACCheelL=ACCheelL(NearleftHSi(j)-HSrefinePre-3:NearleftHSi(j)+HSrefinePost-3);   
    if HSrefinePre>10
        for k=1:length(CutACCheelL)
            if CutACCheelL(k)>=mean(ACCheelL((NearleftHSi(j)-8):(NearleftHSi(j)+2)))+5   
                CutACCheelL(k)=mean(ACCheelL((NearleftHSi(j)-8):(NearleftHSi(j)+2)))+5;
            end
        end
    end
    
    if length(CutACCheelL)>=3    
        [C,D] = findpeaks(CutACCheelL(1:end),'Threshold',0.0001); 
        if length(D)>=1               
            leftHSi(j) = D(1)+NearleftHSi(j)-(HSrefinePre+1);   
        else
            leftHSi(j) = NearleftHSi(j);   
        end
    else
        leftHSi(j) = NearleftHSi(j);   
    end
        
end

% right
[pks,rightHSpeak] = findpeaks(RHEEyF(1:end),'minpeakdistance',150); 
j=0;
for i = 1:length(rightHSpeak)-1   
    in = rightHSpeak(i);
    fi = rightHSpeak(i+1);        
    [A,B] = min(RHEEyF(in:fi));  
    j=j+1;
    NearrightHSi(j) = B+in-1;     
    
    CutACCheelR=ACCheelR(NearrightHSi(j)-HSrefinePre-3:NearrightHSi(j)+HSrefinePost-3);   
    if HSrefinePre>10
        for k=1:length(CutACCheelR)
            if CutACCheelR(k)>=mean(ACCheelR((NearrightHSi(j)-8):(NearrightHSi(j)+2)))+5
                CutACCheelR(k)=mean(ACCheelR((NearrightHSi(j)-8):(NearrightHSi(j)+2)))+5;
            end
        end
    end
    
    if length(CutACCheelR)>=3    
        [C,D] = findpeaks(CutACCheelR(1:end),'Threshold',0.0001);
        if length(D)>=1               
            rightHSi(j) = D(1)+NearrightHSi(j)-(HSrefinePre+1);   
        else
            rightHSi(j) = NearrightHSi(j);   
        end
    else
        rightHSi(j) = NearrightHSi(j);   
    end
end

%% UPDATED PART
% if size(rightHSi,2)>size(leftHSi,2)
%     rightHSi(:,1) = [];    
% end  
%% toe off  
for i = 4:m-3        % m = size(RTOE,1);
    VtoeR(i-3) = (-RTOE(i+2) + 8*RTOE(i+1) - 8*RTOE(i-1) + RTOE(i-2))/(12*h);   
    VtoeL(i-3) = (-LTOE(i+2) + 8*LTOE(i+1) - 8*LTOE(i-1) + LTOE(i-2))/(12*h);
    ACCtoeR(i-3) = (-RTOE(i+2) + 16*RTOE(i+1) - 30*RTOE(i) + 16*RTOE(i-1) - RTOE(i-2) )/(12*h^2);
    ACCtoeL(i-3) = (-LTOE(i+2) + 16*LTOE(i+1) - 30*LTOE(i) + 16*LTOE(i-1) - LTOE(i-2) )/(12*h^2);    
end

% left
j=0;
for i = 1:length(leftHSi)-1   
    in = leftHSi(i)+75;     
    fi = leftHSi(i+1);        
    [pks,locs] = findpeaks(ACCtoeL(in-3:fi-3),'minpeakdistance',TOminpeakdistance,'minpeakheight',TOminpeakheight);  
    j=j+1;
    leftTOi(j) = locs(1)+in-1;   
end

% right
j=0;
for i = 1:length(rightHSi)-1   
    in = rightHSi(i)+75;  
    fi = rightHSi(i+1);       
    [pks,locs] = findpeaks(ACCtoeR(in-3:fi-3),'minpeakdistance',TOminpeakdistance,'minpeakheight',TOminpeakheight);  
    j=j+1;
    rightTOi(j) = locs(1)+in-1; 
end

%% Trim Stride
RHS=rightHSi; LTO=leftTOi; LHS=leftHSi; RTO=rightTOi;
i=0;   

while (length(RHS)~=length(LHS) | length(RTO)~=length(LTO) | length(RHS)~=length(RTO)) & i<20
    if RHS(1)>LTO(1) 
        LTO=LTO(2:end);
    end
    if LTO(1)>LHS(1)
        LHS=LHS(2:end);
    end
    if LHS(1)>RTO(1)
        RTO=RTO(2:end);
    end

    if RHS(end)>LTO(end)
        RHS=RHS(1:end-1);
    end
    if LTO(end)>LHS(end)
        LTO=LTO(1:end-1);
    end
    if LHS(end)>RTO(end)
        LHS=LHS(1:end-1);
    end
    
  
    if RHS(2)<LTO(1)
        RHS=RHS(2:end);
    end
    if LTO(2)<LHS(1)
        LTO=LTO(2:end);
    end
    if LHS(2)<RTO(1)
        LHS=LHS(2:end);
    end

    if RHS(end)<LTO(end-1)
        LTO=LTO(1:end-1);
    end
    if LTO(end)<LHS(end-1)
        LHS=LHS(1:end-1);
    end
    if LHS(end)<RTO(end-1)
        RTO=RTO(1:end-1);
    end
    i=i+1;
end

if length(RHS)~=length(LHS) | length(RTO)~=length(LTO) | length(RHS)~=length(RTO)   
    GE=[];
    GEInMiddleDeleted=[];
    fprintf('\nWarning: Gait events not equal.\n');
    return
end

%% CHECK 1 

% check every row
check=[];
j=0;
for i=1:length(RHS)
    if RHS(i)>LTO(i) | LTO(i)>LHS(i) | LHS(i)>RTO(i)
        j=j+1;
        check(j)=i;
    end
end

X=1;
GEInMiddleDeleted=[];
while length(check)>=1 
    k=check(1);
    if RHS(k)>LTO(k)
        GEInMiddleDeleted(X)=LTO(k); LTO(k)=[]; RHS=RHS(1:end-1); LHS=LHS(1:end-1); RTO=RTO(1:end-1); X=X+1;
    elseif LTO(k)>LHS(k)
        GEInMiddleDeleted(X)=LHS(k); LHS(k)=[]; RHS=RHS(1:end-1); LTO=LTO(1:end-1); RTO=RTO(1:end-1); X=X+1;
    elseif LHS(k)>RTO(k)
        GEInMiddleDeleted(X)=RTO(k); RTO(k)=[]; RHS=RHS(1:end-1); LTO=LTO(1:end-1); LHS=LHS(1:end-1); X=X+1;
    end
    
    check=[];
    j=0;
    for i=1:length(RHS)
        if RHS(i)>LTO(i) | LTO(i)>LHS(i) | LHS(i)>RTO(i)
            j=j+1;
            check(j)=i;
        end
    end
    
end

%% CHECK 2

avemean_RHS = mean(diff(RHS));
X=1;
for i = 1:length(RHS)-1
    if RHS(i+1)-RHS(i) < avemean_RHS - (avemean_RHS/5)
       deleted_RHS(X) = i+1; X=X+1;
    end
end 

if exist('deleted_RHS','var') == 1
RHS(deleted_RHS)=[];LHS(deleted_RHS)=[];RTO(deleted_RHS)=[];LTO(deleted_RHS)=[];
end

avemean_LHS = mean(diff(LHS));
X=1;
for i = 1:length(LHS)-1
    if LHS(i+1)-LHS(i) < avemean_LHS - (avemean_LHS/5)
       deleted_LHS(X) = i+1; X=X+1;
    end
end

if exist('deleted_LHS','var') == 1
RHS(deleted_LHS)=[];LHS(deleted_LHS)=[];RTO(deleted_LHS)=[];LTO(deleted_LHS)=[];
end

avemean_RTO = mean(diff(RTO));
X=1;
for i = 1:length(RTO)-1
    if RTO(i+1)-RTO(i) < avemean_RTO - (avemean_RTO/5)
       deleted_RTO(X) = i+1; X=X+1;
    end
end 

if exist('deleted_RTO','var') == 1
RHS(deleted_RTO)=[];LHS(deleted_RTO)=[];RTO(deleted_RTO)=[];LTO(deleted_RTO)=[];
end 

avemean_LTO = mean(diff(LTO));
X=1;
for i = 1:length(LTO)-1
    if LTO(i+1)-LTO(i) < avemean_LTO - (avemean_LTO/5)
       deleted_LTO(X) = i+1; X=X+1;
    end
end

if exist('deleted_LTO','var') == 1
RHS(deleted_LTO)=[];LHS(deleted_LTO)=[];RTO(deleted_LTO)=[];LTO(deleted_LTO)=[];
end   

    %output GE
    GE=[];
    GE(:,1)=RHS;
    GE(:,2)=LTO;
    GE(:,3)=LHS;
    GE(:,4)=RTO;

% if backwards walking
if nargin == 8 && BW == 1
    GE = flipud(GE);
end
