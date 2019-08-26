function[GE]=get_gaitEvents_GRF(LeftBeltForce,RightBeltForce)
% Inputs are force data from each belt, such as 'FP1For' which contains 3 column vectors
% forces x,y,and z match the original treadmill axes
% +x to the right, +y is vertical, and +z is to the back
% Output is a matrix showes the gait events' index (corresponding to row number of txt file from Dflow)
% the order is:
% RHS LTO LHS RTO
% RHS ...


%%
    % left foot divide
    IS=struct;
    IS.threshold=200;
    IS.Lraw=find(((LeftBeltForce(1:end-1,2)-IS.threshold).*(LeftBeltForce(2:end,2)-IS.threshold))<0);   
        if (LeftBeltForce(IS.Lraw(1)+1,2)-LeftBeltForce(IS.Lraw(1),2))>0    
            IS.L=IS.Lraw(2:end);
        else
            IS.L=IS.Lraw;
        end
        if mod(length(IS.L),2)==1  
            IS.L=IS.L(1:end-1);
        end
    IS.Ls=IS.L(1:2:end);   
    IS.Le=IS.L(2:2:end);   

    %right foot divide
    IS.Rraw=find(((RightBeltForce(1:end-1,2)-IS.threshold).*(RightBeltForce(2:end,2)-IS.threshold))<0);   
        if (RightBeltForce(IS.Rraw(1)+1,2)-RightBeltForce(IS.Rraw(1),2))>0    
            IS.R=IS.Rraw(2:end);
        else
            IS.R=IS.Rraw;
        end
        if mod(length(IS.R),2)==1  
            IS.R=IS.R(1:end-1);
        end
    IS.Rs=IS.R(1:2:end);   
    IS.Re=IS.R(2:2:end);   

    % find toe off and heel strike for left foot
    GEraw=struct;
    GEraw.thresholdfix=10;  
    for i=1:length(IS.Ls)   
        GEraw.Lthresholdfloat(i)=min(LeftBeltForce(IS.Ls(i):IS.Le(i),2))+10;  
    end
    GEraw.Lthresholdfloat=GEraw.Lthresholdfloat';
    GEraw.LTO=[];
    GEraw.LHS=[];
    for i=1:length(IS.Ls)
        %GEraw.Lbelowthreshold=find(LeftBeltForce(IS.Ls(i):IS.Le(i),2)<GEraw.thresholdfix);        
        GEraw.Lbelowthreshold=find(LeftBeltForce(IS.Ls(i):IS.Le(i),2)<GEraw.Lthresholdfloat(i));   
        GEraw.LTO=union(IS.Ls(i)+GEraw.Lbelowthreshold(1)-1,GEraw.LTO);
        GEraw.LHS=union(IS.Ls(i)+GEraw.Lbelowthreshold(end)-1,GEraw.LHS);
    end
    GEraw.LTO=GEraw.LTO';
    GEraw.LHS=GEraw.LHS';
    
    % find toe off and heel strike for right foot
    for i=1:length(IS.Rs)   
        GEraw.Rthresholdfloat(i)=min(RightBeltForce(IS.Rs(i):IS.Re(i),2))+10;  
    end
    GEraw.Rthresholdfloat=GEraw.Rthresholdfloat';
    GEraw.RTO=[];
    GEraw.RHS=[];
    for i=1:length(IS.Rs)
        %GEraw.Rbelowthreshold=find(RightBeltForce(IS.Rs(i):IS.Re(i),2)<GEraw.thresholdfix);        
        GEraw.Rbelowthreshold=find(RightBeltForce(IS.Rs(i):IS.Re(i),2)<GEraw.Rthresholdfloat(i));   
        GEraw.RTO=union(IS.Rs(i)+GEraw.Rbelowthreshold(1)-1,GEraw.RTO);
        GEraw.RHS=union(IS.Rs(i)+GEraw.Rbelowthreshold(end)-1,GEraw.RHS);
    end
    GEraw.RTO=GEraw.RTO';
    GEraw.RHS=GEraw.RHS';
    
    % arrange GE in order
    GEraw.GE1=union(GEraw.LTO,GEraw.LHS); 
    GEraw.GE1=union(GEraw.GE1,GEraw.RTO);
    GEraw.GE1=union(GEraw.GE1,GEraw.RHS);
    GEraw.GE2=GEraw.GE1;
    for i=1:3     
        if GEraw.GE2(1)~=GEraw.RHS(1)
            GEraw.GE2=GEraw.GE2(2:end);
        end
    end
    GEraw.GE3=GEraw.GE2;
    for i=1:3     
        if GEraw.GE3(end)~=GEraw.RTO(end)
            GEraw.GE3=GEraw.GE3(1:end-1);
        end
    end
%     GE(:,1)=GEraw.GE3(1:2:end);
%     GE(:,2)=GEraw.GE3(2:2:end);

if mod(numel(GEraw.GE3),4) ~= 0
    GEraw.GE3 = GEraw.GE3(1:end-mod(numel(GEraw.GE3),4));
end

GE1 = reshape(GEraw.GE3,[4, numel(GEraw.GE3)/4]);
GE = GE1';