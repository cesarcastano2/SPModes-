% close all
clearvars
clc
% % 
% subjs = {'Ll01' 'Ll02' 'Ll03' 'Ll04' 'Ll05' 'Ll06' 'Ll07' 'Ll10'};
subjs = {'Ll03'};
% conds = {'level_050' 'level_075' 'level_100' 'level_125' 'level_selfpaced' 'incline_selfpaced' 'decline_selfpaced' 'incline_075' 'decline_075'};
conds = {'level_selfpaced'};
theta = [-9 0 9];                                                       
projfolder = pwd;

color = [0 0 1; 0 1 0; 1 0 0];
pitchval = [0 9 -9];
for s = 1:length(subjs)
    
    for c = 1:length(conds)
        %% LOAD DATA FILE
        
        % if you want, you can set the file to load here instead of using the gui
        dflow_file = [projfolder '/' subjs{s} '/' subjs{s} '_' conds{c} '0001.txt'];
        dflow_treadmill_file = [projfolder '/' subjs{s} '/' subjs{s} '_' conds{c} '_treadmill0001.txt'];
        
        if ~exist('dflow_file','var')
            [FileName,PathName,FilterIndex]  = uigetfile('*.txt');
            dflow_file = [PathName FileName];
        end
        
        if ~exist('dflow_treadmill_file','var')
            [FileName,PathName,FilterIndex]  = uigetfile('*treadmill0001.txt');
            dflow_treadmill_file = [PathName FileName];
        end
        
        % load DFLOW data (markers & forces from treadmill)
        % FP1 = left, FP2 = right
        [Frame_df, Time_df, markers_df, forces_df, startidx, stopidx, Total] = import_dflow(dflow_file);
        
        datatreadmill = import_dflow_treadmill(dflow_treadmill_file);
%         Time_treadmill = datatreadmill.Time-datatreadmill.Time(1);
        
        full_data_treadmill = datatreadmill; 
    
       %% Calc Speed from treadmill

datatreadmill.Time = interpft(datatreadmill.Time(:,:),length(Frame_df));
datatreadmill.LeftBeltSpeed = interpft(datatreadmill.LeftBeltSpeed(:,:),length(Frame_df));
datatreadmill.RightBeltSpeed  = interpft(datatreadmill.RightBeltSpeed(:,:),length(Frame_df));
datatreadmill.Pitch = interpft(datatreadmill.Pitch(:,:),length(Frame_df));

datatreadmill.Time = datatreadmill.Time(startidx:stopidx);
datatreadmill.LeftBeltSpeed =  datatreadmill.LeftBeltSpeed(startidx:stopidx);
datatreadmill.RightBeltSpeed  = datatreadmill.RightBeltSpeed(startidx:stopidx);
datatreadmill.Pitch = datatreadmill.Pitch(startidx:stopidx);

  
        markers_df.LASI = markers_df.LASI(startidx:stopidx,:);
        markers_df.RASI = markers_df.RASI(startidx:stopidx,:);
        markers_df.LPSI = markers_df.LPSI(startidx:stopidx,:);
        markers_df.RPSI = markers_df.RPSI(startidx:stopidx,:);
        markers_df.LKNE = markers_df.LKNE(startidx:stopidx,:);
        markers_df.LTHI = markers_df.LTHI(startidx:stopidx,:);
        markers_df.LANK = markers_df.LANK(startidx:stopidx,:);
        markers_df.LTIB = markers_df.LTIB(startidx:stopidx,:);
        markers_df.LTOE = markers_df.LTOE(startidx:stopidx,:);
        markers_df.LHEE = markers_df.LHEE(startidx:stopidx,:);
        markers_df.RKNE = markers_df.RKNE(startidx:stopidx,:);
        markers_df.RTHI = markers_df.RTHI(startidx:stopidx,:);
        markers_df.RANK = markers_df.RANK(startidx:stopidx,:);
        markers_df.RTIB = markers_df.RTIB(startidx:stopidx,:);
        markers_df.RTOE = markers_df.RTOE(startidx:stopidx,:);
        markers_df.RHEE = markers_df.RHEE(startidx:stopidx,:);
        
        forces_df.FP1Cop = forces_df.FP1Cop(startidx:stopidx,:);
        forces_df.FP1For = forces_df.FP1For(startidx:stopidx,:);
        forces_df.FP1Mom = forces_df.FP1Mom(startidx:stopidx,:);
        forces_df.FP2Cop = forces_df.FP2Cop(startidx:stopidx,:);
        forces_df.FP2For = forces_df.FP2For(startidx:stopidx,:);
        forces_df.FP2Mom = forces_df.FP2Mom(startidx:stopidx,:);
        
        Frame_df = Frame_df(startidx:stopidx,:);
        Time_df = Time_df(startidx:stopidx,:);
        


         
%% FINDING START 
start_names = {'Ll01level_selfpaced','Ll01incline_selfpaced', 'Ll01decline_selfpaced', 'Ll02level_selfpaced', 'Ll02incline_selfpaced', 'Ll02decline_selfpaced', 'Ll03level_selfpaced', 'Ll03incline_selfpaced', 'Ll03decline_selfpaced', 'Ll04level_selfpaced', 'Ll04incline_selfpaced', 'Ll04decline_selfpaced', 'Ll05level_selfpaced', 'Ll05incline_selfpaced', 'Ll05decline_selfpaced', 'Ll06level_selfpaced', 'Ll06incline_selfpaced', 'Ll06decline_selfpaced', 'Ll07level_selfpaced', 'Ll07incline_selfpaced', 'Ll07decline_selfpaced', 'Ll10level_selfpaced', 'Ll10incline_selfpaced', 'Ll10decline_selfpaced'};
start_values = [327 328 327 283 308 289 292 292 291 279 290 288 294 296 295 292 294 296 287 286 286 296 333 297];
start_location = find(strcmp([subjs{s} conds{c}],start_names));
start_value_pref = start_values(start_location);

        %% HARD START / STOP
        Dflow_size = length(Time_df);
%         Start_Time = find(Time_df>5);
%         Start_Time = Start_Time(1);

        Z=1;
%         
        if Z==sum(strcmp(subjs(s),'Ll08')) && Z==sum(strcmp(conds(c),'decline_selfpaced')) 
           Stop_Time_min = find(Time_df > floor(Time_df(end)));
           Stop_Time = Stop_Time(1);
        else
            
           Stop_Time_min = floor(Time_df(end));
           Stop_Time = find(Time_df > Stop_Time_min);
           Stop_Time = Stop_Time(1);
        end
        Stop_Time = Stop_Time(end);
        
%         Start_Time = find(Time_df >(Stop_Time_min - 60));
        Start_Time = find(Time_df >(Stop_Time_min - start_value_pref));
        Start_Time = Start_Time(1);
      
%          Stop_Time = find(Time_df >(Stop_Time_min - start_value_pref)+275);
%          Stop_Time = Stop_Time(1);


%         if Z==sum(strcmp(subjs(s),'Ll08')) && Z==sum(strcmp(conds(c),'decline_selfpaced')) 
%            Stop_Time = find(Time_df<290);
%         else
%            Stop_Time = find(Time_df<300);
%         end
%         Stop_Time = Stop_Time(end);
        
%         for i = length(markers_df.labels)
%        eval(['markers_df.' markers_df.labels(i) ' = (Start_Time(startidx:stopidx,markers.idx(m:m+2)));']); 
        
        markers_df.LASI = markers_df.LASI(Start_Time:Stop_Time,:);
        markers_df.RASI = markers_df.RASI(Start_Time:Stop_Time,:);
        markers_df.LPSI = markers_df.LPSI(Start_Time:Stop_Time,:);
        markers_df.RPSI = markers_df.RPSI(Start_Time:Stop_Time,:);
        markers_df.LKNE = markers_df.LKNE(Start_Time:Stop_Time,:);
        markers_df.LTHI = markers_df.LTHI(Start_Time:Stop_Time,:);
        markers_df.LANK = markers_df.LANK(Start_Time:Stop_Time,:);
        markers_df.LTIB = markers_df.LTIB(Start_Time:Stop_Time,:);
        markers_df.LTOE = markers_df.LTOE(Start_Time:Stop_Time,:);
        markers_df.LHEE = markers_df.LHEE(Start_Time:Stop_Time,:);
        markers_df.RKNE = markers_df.RKNE(Start_Time:Stop_Time,:);
        markers_df.RTHI = markers_df.RTHI(Start_Time:Stop_Time,:);
        markers_df.RANK = markers_df.RANK(Start_Time:Stop_Time,:);
        markers_df.RTIB = markers_df.RTIB(Start_Time:Stop_Time,:);
        markers_df.RTOE = markers_df.RTOE(Start_Time:Stop_Time,:);
        markers_df.RHEE = markers_df.RHEE(Start_Time:Stop_Time,:);
        
        forces_df.FP1Cop = forces_df.FP1Cop(Start_Time:Stop_Time,:);
        forces_df.FP1For = forces_df.FP1For(Start_Time:Stop_Time,:);
        forces_df.FP1Mom = forces_df.FP1Mom(Start_Time:Stop_Time,:);
        forces_df.FP2Cop = forces_df.FP2Cop(Start_Time:Stop_Time,:);
        forces_df.FP2For = forces_df.FP2For(Start_Time:Stop_Time,:);
        forces_df.FP2Mom = forces_df.FP2Mom(Start_Time:Stop_Time,:);
        
        Frame_df = Frame_df(Start_Time:Stop_Time,:);
        Time_df = Time_df(Start_Time:Stop_Time,:);
        
        datatreadmill.Time = datatreadmill.Time(Start_Time:Stop_Time,:);
        datatreadmill.LeftBeltSpeed =  datatreadmill.LeftBeltSpeed(Start_Time:Stop_Time,:);
        datatreadmill.RightBeltSpeed  = datatreadmill.RightBeltSpeed(Start_Time:Stop_Time,:); 
        datatreadmill.Pitch = datatreadmill.Pitch(Start_Time:Stop_Time,:);
        
        conds_all = {'level_050' 'level_075'  'level_100' 'level_125' 'level_selfpaced' 'decline_selfpaced' 'decline_075' 'incline_selfpaced' 'incline_075'};
        
      
        %%
        % ML, horizontal/parallel, vertical/perpendicular
        %         F_parallel = forces_df_c.FP2For(:,2)*cosd(theta(c)) + forces_df_c.FP2For(:,3)*sind(theta(c));
        %         F_perpendicular = forces_df_c.FP2For(:,2)*sind(theta(c)) + forces_df_c.FP2For(:,3)*cosd(theta(c));
        %         Rgrf = [forces_df_c.FP2For(:,1) F_parallel F_perpendicular];
        %
        %         F_parallel = forces_df_c.FP1For(:,2)*cosd(theta(c)) + forces_df_c.FP1For(:,3)*sind(theta(c));
        %         F_perpendicular = forces_df_c.FP1For(:,2)*sind(theta(c)) + forces_df_c.FP1For(:,3)*cosd(theta(c));
        %         Lgrf = [forces_df_c.FP1For(:,1) F_parallel F_perpendicular];
        
      
        %% GET GAIT EVENTS
        HSrefinePre=10;
        HSrefinePost= 5;
        TOminpeakheight=4 ;
        TOminpeakdistance=40;
        BW=0;
        Time = Time_df;
        markers4GE = {'RHEE' 'LHEE' 'RANK' 'LANK' 'RTOE' 'LTOE'};
%         markers_df.labels(1,:) = {'LASIS', 'RASIS', 'LPSIS', 'RPSIS', 'LKNE', 'LTHI', 'LANK', 'LTIB', 'LTOE', 'LHEE', 'RKNE', 'RTHI', 'RANK', 'RTIB', 'RTOE', 'RHEE'};
        llmarkers= markers_df;
        Z=1;
        if Z==sum(strcmp(conds(c),conds_all)) 
        [RHS,LTO,LHS,RTO,GE,GEInMiddleDeleted] = GaitEvents_allslopes(Time,llmarkers,markers4GE, HSrefinePre, HSrefinePost, TOminpeakdistance, TOminpeakheight,BW); 
        end 
        
        if Z==sum(strcmp(conds(c),conds_all)) 
        [RHS,LTO,LHS,RTO,GE,GEInMiddleDeleted] = GaitEvents_allslopes(Time,llmarkers,markers4GE, HSrefinePre, HSrefinePost, TOminpeakdistance, TOminpeakheight,BW); 
        end 
        
        RHS = RHS';
        LTO = LTO';
        LHS = LHS';
        RTO = RTO';
        GEgood=GE;
       
%         acceptable = GEgood(:,1)>7200;
%         GEgood = GEgood(acceptable(:),:);
        PPP=[];
        leftheel_GE = [];
        
        PPP=repmat(GEgood(2:end,1),1);
        leftheel_GE = repmat(GEgood(2:end,3),1);
        GEgood(end,:)=[];
        GEgood(:,5)=PPP;
      
         
%% Calc Speed from treadmill


blockSize = [240, 1];
meanFilterFunction = @(theBlockStructure) mean2(theBlockStructure.data(:));
blockAveragedDownSignal = blockproc(datatreadmill.RightBeltSpeed(:,1), blockSize, meanFilterFunction);
[rows, columns] = size(blockAveragedDownSignal);


yy = find(full_data_treadmill.RightBeltSpeed(1:30000)==0,1,'last');
% yy = yy(1);

% x=(1:length(full_data_treadmill.RightBeltSpeed(1:Stop_Time,1)));
% y=full_data_treadmill.RightBeltSpeed(1:Stop_Time,1);

x=(1:length(full_data_treadmill.RightBeltSpeed(yy:(yy+75000),1)));
y=full_data_treadmill.RightBeltSpeed(yy:yy+75000,1);
Speed_Time = full_data_treadmill.Time(yy:(yy+75000),1);
% x = full_data_treadmill.Time(yy:(yy+75000),1);
Speed_Time = Speed_Time - Speed_Time(1);
x = Speed_Time - Speed_Time(1);


Ll01_speed(:,1) = Speed_Time;
Ll01_speed(:,c+1) = y;
% 

fitpoints = find(y==0,1,'last');

% fitpoints = fitpoints(5000);
% 
% w = x(fitpoints:end); 
% z = y(fitpoints:end);

% figure(s)
% hold on
% scatter(x,y,'.','DisplayName',[subjs{s} ' ' conds{c}],'MarkerEdgeColor',color(c,:))
% legend()
x1 = [];
x1 = x(fitpoints:end)- x(fitpoints);
y1 = [];
y1 = y(fitpoints:end);

if c == 1 
full_treadmill_speed.level(s,:) = y1;
end
if c == 2 
full_treadmill_speed.incline(s,:) = y1;
end
if c == 3
full_treadmill_speed.decline(s,:) = y1;
end



 %% Convert to conventional coords
        for m = 1:length(markers_df.labels)
            eval(['markers_df_c.' markers_df.labels{m} ' = convert_coords2conventional(markers_df. ' markers_df.labels{m} ');']);
        end
        
        for m = 1:length(forces_df.labels)
            eval(['forces_df_c.' forces_df.labels{m} ' = convert_coords2conventional(forces_df. ' forces_df.labels{m} ');']);
        end
        
%% step length / Step Time / Stride Time 
stride_speed = [];
stride_time = [];
steplength_speed = [];
steplength_time = [];
stridelength_value = [];
stridewidth_value = [];


        for i = 1:length(GEgood(:,1))
            steplength_speed(i,1) = nanmean(datatreadmill.RightBeltSpeed(GEgood(i,1):GEgood(i,3)));
            steplength_speed(i,2) = nanmean(datatreadmill.LeftBeltSpeed(GEgood(i,3):GEgood(i,5)));
        end  

        for i = 1:length(GEgood(:,1))
            steplength_time(i,1) = range(Time_df(GEgood(i,1):GEgood(i,3)));
            steplength_time(i,2) = range(Time_df(GEgood(i,3):GEgood(i,5)));
        end
        
        for i = 1:length(GEgood(:,1))
            stride_time(i,1) = range(Time_df(GEgood(i,1):GEgood(i,5)));
            stride_time(i,2) = range(Time_df(GEgood(i,3):leftheel_GE(i,1)));
        end
        
      
        for i = 1:length(GEgood(:,1))
            stride_speed(i,1) = nanmean(datatreadmill.RightBeltSpeed(GEgood(i,1):GEgood(i,5)));
            stride_speed(i,2) = nanmean(datatreadmill.LeftBeltSpeed(GEgood(i,3):leftheel_GE(i,1)));    
        end
        
        stridelength_value(:,1) = stridelength(markers_df_c.RHEE(:,2), markers_df_c.RHEE(:,2), GEgood(:,5), GEgood(:,1),stride_speed(:,1), stride_time(:,1),pitchval(c));
        stridelength_value(:,2) = stridelength(markers_df_c.LHEE(:,2), markers_df_c.LHEE(:,2), leftheel_GE(:,1), GEgood(:,3),stride_speed(:,2), stride_time(:,2),pitchval(c));
        
        stridewidth_value(:,1) = stridewidth(markers_df_c.RHEE(:,1), markers_df_c.LHEE(:,1), GEgood(:,1), GEgood(:,3));
        stridewidth_value(:,2) = abs(stridewidth(markers_df_c.LHEE(:,1),markers_df_c.RHEE(:,1),  GEgood(:,3), GEgood(:,5)));
        %%
%         h(c)=figure;
%         GEgood = gaiteventCheck3(RHS, LTO, LHS, RTO, forces_df_c.FP2For(:,3), forces_df_c.FP1For(:,3));

        figure
        
        plot(forces_df_c.FP2For(:,3),'r');
        hold on
        plot(forces_df_c.FP1For(:,3),'b');
        %         plot(min(forces_df_c.FP2For(:,3), forces_df_c.FP1For(:,3)), 'k', 'linewidth', 2)
        %         plot(min(forces_df_c.FP2For(:,3), forces_df_c.FP1For(:,3)), 'k')
        plot(GEgood(:,1), zeros(size(GEgood(:,1))), 'rx', GEgood(:,2), zeros(size(GEgood(:,1))), 'bo', GEgood(:,3), zeros(size(GEgood(:,1))), 'bx', GEgood(:,4), zeros(size(GEgood(:,1))), 'ro', GEgood(:,5), zeros(size(GEgood(:,1))), 'rx')
        ylimits = ylim(gca);
        title([subjs{s} ', ' conds{c}])
        %             axis([0 RHS(10) ylimits(1) ylimits(2)]);
        
        
%         savefig(h,'Incline_gait_ACC.fig')
%        close(gcf) 
        %% SAVE GE VARIABLES
        
%         filename = [subjs(s) conds(c)];
%         filename = strjoin(filename,{'_'});
%         filename =  strcat(filename,'.mat');
%         save(filename,'GE')
        
        %% Step width & Step length
        
        clear sl sw badidx1 badidx2 allbadsl
        
%         sl(:,1) = steplength(markers_df_c.RHEE(:,2), markers_df_c.LHEE(:,2), GEgood(:,1));
%         sw(:,1) = stepwidth(markers_df_c.RHEE(:,1), markers_df_c.LHEE(:,1), GEgood(:,1));
%         
%         sl(:,2) = steplength(markers_df_c.LHEE(:,2), markers_df_c.RHEE(:,2), GEgood(:,3));
%         sw(:,2) = stepwidth(markers_df_c.LHEE(:,1), markers_df_c.RHEE(:,1), GEgood(:,3)); 
        sl=[]; sw=[]; sp=[];  st=[];
        
        sl(:,1) = steplength(markers_df_c.LHEE(:,2), markers_df_c.RHEE(:,2), GEgood(:,3), GEgood(:,1),steplength_speed(:,1), steplength_time(:,1),pitchval(c));
        sw(:,1) = stepwidth(markers_df_c.RHEE(:,1), markers_df_c.LHEE(:,1), GEgood(:,3));
        
        sl(:,2) = steplength(markers_df_c.RHEE(:,2), markers_df_c.LHEE(:,2), GEgood(:,5), GEgood(:,3),steplength_speed(:,2), steplength_time(:,2),pitchval(c));
        sw(:,2) = stepwidth(markers_df_c.LHEE(:,1), markers_df_c.RHEE(:,1), GEgood(:,5));
        goodidx1 = find(sl(:,1) > 0.9*mean(sl(:,1)));
        goodidx2 = find(sl(:,2) > 0.9*mean(sl(:,2)));
        goodsl = intersect(goodidx1, goodidx2); 
        tempsl = sl(goodsl,:);
        tempsw = sw(goodsl,:);
        stepcount = length(sl);
%         sl = tempsl;
%         sw = tempsw;
        sp = datatreadmill.LeftBeltSpeed(:);
%         sp = avespeed(39:end)';
        for sn = 1:length(GEgood)
            st(sn,1) = Time_df(GEgood(sn,5))-Time_df(GEgood(sn,1));
        end
%         for sn = 41:length(GEgood)
%             st(sn,1) = Time_df(GEgood(sn,5))-Time_df(GEgood(sn,1));
%         end
%         
        
       
        
        sp_length(s,c) = length(sp);
        sp_plot(s,c) = mean(sp);
        
        
%         sp_length(s,c) = length(sp)/2;
%         sp_1half(s,c) = mean(sp(1:sp_length(s,c)));
%         sp_2half(s,c) = mean(sp(sp_length(s,c)+1:end));
        
        spvar_length(s,c) = length(sp)/2;
        spvar_1half(s,c) = std(sp(1:spvar_length(s,c)));
        spvar_2half(s,c) = std(sp(spvar_length(s,c)+1:end));
        
        
        st_length(s,c) = length(st)/2;
        st_1half(s,c) = std(st(1:st_length(s,c)));
        st_2half(s,c) = std(st(st_length(s,c)+1:end));
        
        stvar_length(s,c) = length(st)/2;
        stvar_1half(s,c) = mean(st(1:stvar_length(s,c)));
        stvar_2half(s,c) = mean(st(stvar_length(s,c)+1:end));
        
        sl_length(s,c) = length(sl(:,1))/2;
        sl_1half(s,c) = mean(sl(1:sl_length(s,c),1));
        sl_2half(s,c) = mean(sl(sl_length(s,c)+1:end,1));
        
        slvar_length(s,c) = length(sl(:,1))/2;
        slvar_1half(s,c) = mean(sl(1:stvar_length(s,c),1));
        slvar_2half(s,c) = mean(sl(stvar_length(s,c)+1:end,1));
       


        GRP.sp(s,c) = mean(sp);
        GRP.st(s,c) = mean(st);
        GRP.sl(s,:,c) = mean(sl);
        GRP.sw(s,:,c) = mean(sw);
        
        GRP.spvar(s,c) = std(sp);
        GRP.stvar(s,c) = std(st);
        GRP.slvar(s,:,c) = std(sl);
        GRP.swvar(s,:,c) = std(sw);



        %% GRF
        
        numstrides = size(GEgood,1);
        npts = 250;
        GEn = 1; %RHS
        
        [time_n, GRF_n, data_nft] = normalize_gaitcycle2(Time_df, forces_df_c.FP2For, GEgood, npts);
        
        temp = nanmean(GRF_n,3);
        GRP.GRF_ML_R(:,s,c) = temp(:,1);
        GRP.GRF_AP_R(:,s,c) = temp(:,2);
        GRP.GRF_v_R(:,s,c) = temp(:,3);
        
        
        %% Calculate segment angles for strides 10-15
        %
        % trunk, thigh, shank, foot
        %
        RTHIminusRKNE = markers_df_c.RTHI - markers_df_c.RKNE;
        RTIBminusRANK = markers_df_c.RTIB - markers_df_c.RANK;
        RTOEminusRHEE = markers_df_c.RTOE - markers_df_c.RHEE;
        RHEEminusRTOE = markers_df_c.RHEE - markers_df_c.RTOE;
        
        angle.RTHI = atan2d(RTHIminusRKNE(:,3),RTHIminusRKNE(:,2));
        angle.RSHA = atan2d(RTIBminusRANK(:,3),RTIBminusRANK(:,2));
        angle.RFOO = atan2d(RTOEminusRHEE(:,3),RTOEminusRHEE(:,2));
        % angle.RFOO = atan2d(RHEEminusRTOE(:,3),RHEEminusRTOE(:,2)); %discont
        
        %%
        
        [time_n, angle_n.RTHI, data_nft] = normalize_gaitcycle2(Time_df, angle.RTHI, GEgood, npts);
        GRP.RTHI(:,s,c) = mean(angle_n.RTHI,2);
        
        [time_n, angle_n.RSHA, data_nft] = normalize_gaitcycle2(Time_df, angle.RSHA, GEgood, npts);
        GRP.RSHA(:,s,c) = mean(angle_n.RSHA,2);
        
        [time_n, angle_n.RFOO, data_nft] = normalize_gaitcycle2(Time_df, angle.RFOO, GEgood, npts);
        GRP.RFOO(:,s,c) = mean(angle_n.RFOO,2);
        
        %% Joint angles
        %
        % knee, ankle
        %
        jangle.RKNE = angle.RTHI-angle.RSHA;
        jangle.RANK = 90-(angle.RSHA-angle.RFOO);
        
        [time_n, jangle_n.RKNE, data_nft] = normalize_gaitcycle2(Time_df, jangle.RKNE, GEgood, npts);
        GRP.RKNE(:,s,c) = nanmean(jangle_n.RKNE,2);
        
        [time_n, jangle_n.RANK, data_nft] = normalize_gaitcycle2(Time_df, jangle.RANK, GEgood, npts);
        GRP.RANK(:,s,c) = nanmean(jangle_n.RANK,2);
        %% PLOTTING JOINT ANGLES 
%         % ankle 
%         for c = 1:length(conds)
%             aveRKNE(:,c) = nanmean([squeeze(GRP.RKNE(:,1,c)) squeeze(GRP.RKNE(:,2,c)) squeeze(GRP.RKNE(:,3,c)) squeeze(GRP.RKNE(:,4,c)) squeeze(GRP.RKNE(:,5,c)) squeeze(GRP.RKNE(:,6,c)) squeeze(GRP.RKNE(:,7,c)) squeeze(GRP.RKNE(:,8,c)) squeeze(GRP.RKNE(:,9,c))],2);
%             aveRANK(:,c) = nanmean([squeeze(GRP.RANK(:,1,c)) squeeze(GRP.RANK(:,2,c)) squeeze(GRP.RANK(:,3,c)) squeeze(GRP.RANK(:,4,c)) squeeze(GRP.RANK(:,5,c)) squeeze(GRP.RANK(:,6,c)) squeeze(GRP.RANK(:,7,c)) squeeze(GRP.RANK(:,8,c)) squeeze(GRP.RANK(:,9,c))],2);
%         end
%         figure()
%         plot(ankle_avg(1,:))
%         
%         % knee
%         for i=1:length(jangle_n.RKNE)
%             knee_avg(i) = nanmean(jangle_n.RKNE(i,:));
%         end
%         figure()
%         plot(ankle_knee(1,:)) 


% figure(1);
% hold on
% plot(st(:,1))
% title('st incline')
% 
% figure(2);
% hold on 
% plot(sl(:,1))
% title('sl incline')
% 
% figure(3);
% hold on 
% plot(sp(:,1))
% title('sp incline')
% 
% figure(4);
% hold on 
% plot(sw(:,1))
% title('sw incline')
% 
% h(1)=figure(1);
% h(2)=figure(2);
% h(3)=figure(3);
% h(4)=figure(4);
% 
% savefig(h,['decline_selfpaced.fig'])

% savefig(h,['level_selfpaced.fig']) 
%% PLOTS VS STEPS
 
sl_all = [];
for i = 1:size(sl,1)
    sl_all = [sl_all sl(i,:)];
end

sw_all = [];
for i = 1:size(sw,1)
    sw_all = [sw_all sw(i,:)];
end

% steplength_time = steplength_time';
steplength_time_all = [];
for i = 1:size(steplength_time,1)
    steplength_time_all = [steplength_time_all steplength_time(i,:)];
end

% steplength_speed = steplength_speed';
steplength_speed_all = [];
for i = 1:size(steplength_speed,1)
    steplength_speed_all = [steplength_speed_all steplength_speed(i,:)];
end

stride_speed_all = [];
for i = 1:size(stride_speed,1)
    stride_speed_all = [stride_speed_all stride_speed(i,:)];
end

stride_time_all = [];
for i = 1:size(stride_time,1)
    stride_time_all = [stride_time_all stride_time(i,:)];
end

stridelength_all = [];
for i = 1:size(stridelength_value,1)
    stridelength_all = [stridelength_all stridelength_value(i,:)];
end

stridewidth_all = [];
for i = 1:size(stridewidth_value,1)
    stridewidth_all = [stridewidth_all stridewidth_value(i,:)];
end


%% PLOT for each subject BEST FIT

    
% v = 0.5:0.1:1.8;
%     figure(100)
%     hold on
%     plot(stride_speed_all,stridelength_all,'.','DisplayName',[subjs{s} ' ' conds{c}],'Color',color(c,:))
% %     legend()
%     xlabel('Speed')
%     ylabel('Stride length')
%     [P,fh] = fitsinglemodelprocess_sl(stridelength_all(1:400),stride_speed_all(1:400));
% % [P,fh] = fitsinglemodelprocess_sl(length_strides.level(s,:),speed_strides.level(s,:));
%     plot(v,fh(v,P),'Color',color(c,:));
%     xlim([0.4 1.8])
%     ylim([0.6 1.8])
%     if c == 1 
%     p.level(s,:) = P;
%     end
%     if c == 2 
%     p.incline(s,:) = P;
%     end
%     if c == 3
%     p.decline(s,:) = P;
%     end
%      
     
%      
%      if c == 1
%      variation_strides.level.slminusfit(s,:) = stridelength_all(1:400) - fh(stride_speed_all(1:400),P);
%      variation_strides.level.avevar(s,:) = mean(var(abs(variation_strides.slminusfit(s,:))));
%      variation_strides.level.stdvar(s,:) = std(var(abs(variation_strides.slminusfit(s,:))));
%      variation_strides.level.totalvar(s,:) = var(stridelength_all(1:400));
%      end
%      if c == 2
%      variation_strides.incline.slminusfit(s,:) = stridelength_all(1:400) - fh(stride_speed_all(1:400),P);
%      variation_strides.incline.avevar(s,:) = mean(var(abs(variation_strides.slminusfit(s,:))));
%      variation_strides.incline.stdvar(s,:) = std(var(abs(variation_strides.slminusfit(s,:))));
%      variation_strides.incline.totalvar(s,:) = var(stridelength_all(1:400));
%      end
%      if c == 3
%      variation_strides.decline.slminusfit(s,:) = stridelength_all(1:400) - fh(stride_speed_all(1:400),P);
%      variation_strides.decline.avevar(s,:) = mean(var(abs(variation_strides.slminusfit(s,:))));
%      variation_strides.decline.stdvar(s,:) = std(var(abs(variation_strides.slminusfit(s,:))));
%      variation_strides.decline.totalvar(s,:) = var(stridelength_all(1:400));
%      end
%      
     
%% PCI METHOD
% [PCI_percent_value,step_stride_norm_value,PCI_CV_OP1_value,PCI_CV_OP2_value,ROI_value]  = PCImethod(steplength_time(1:100,:),stride_time(1:100,:));
% 
% if c == 1
%   PCI_values.level.steadystatevalue(s) = ROI_value;
% end
% if c == 2      
%   PCI_values.incline.steadystatevalue(s) = ROI_value;      
% end
% if c == 3
%   PCI_values.decline.steadystatevalue(s) = ROI_value;        
% end
% PCI_steadystate_value = [];

%% STD Method 
[STV_OP1_value,STV_OP1_median_value,STV_OP1_steadypoint_value_all,STV_OP1_steadypoint_value] = stdmethod(stride_speed_all(1:400));

if c == 1
  STV_OP1.level.steadystatevalue(:,s) = STV_OP1_steadypoint_value;
  STV_OP1.level.steadystatevalue_all(:,s) = length(STV_OP1_steadypoint_value_all)/400;
  STV_OP1.level.medianvalue(:,s) = STV_OP1_median_value;
  STV_OP1.level.STDvalues(:,s) = STV_OP1_value;
end
if c == 2      
  STV_OP1.incline.steadystatevalue(:,s) = STV_OP1_steadypoint_value;
  STV_OP1.incline.steadystatevalue_all(:,s) = length(STV_OP1_steadypoint_value_all)/400;
  STV_OP1.incline.medianvalue(:,s) = STV_OP1_median_value;
  STV_OP1.incline.STDvalues(:,s) = STV_OP1_value;
end
if c == 3
  STV_OP1.decline.steadystatevalue(:,s) = STV_OP1_steadypoint_value; 
  STV_OP1.decline.steadystatevalue_all(:,s) = length(STV_OP1_steadypoint_value_all)/400;
  STV_OP1.decline.medianvalue(:,s) = STV_OP1_median_value;
  STV_OP1.decline.STDvalues(:,s) = STV_OP1_value;
end
% STV_OP1_steadypoint_value = [];
% STV_OP1_value = [];
% STV_OP1_median_value = [];
% STV_OP1_steadypoint_value_all = [];

%% STD method plots
% 
% PLOT FOR SINGLE SUBJECT
figure 
title([subjs{s} ' ' conds{c}])
yyaxis left
plot(1:length(stride_speed_all(1:394)),stride_speed_all(1:394),'.-'); 
xlim([0 394])
hold on 
scatter(STV_OP1_steadypoint_value_all,stride_speed_all(STV_OP1_steadypoint_value_all),'.r')
ylabel('Stride Speed')
yyaxis right 
plot(STV_OP1.level.STDvalues((1:394),s),'.-');
ylim([0 0.15])
ylabel('STD Values')
xlabel('Stride Number')
hold on
% yline(STV_OP1.level.medianvalue(:,s),'r');
% xline(STV_OP1.level.steadystatevalue(:,s));
patch('Faces',[1 2 3 4],'Vertices',[0 0;0 STV_OP1.level.medianvalue(:,s);394 STV_OP1.level.medianvalue(:,s);394 0],'FaceColor','red','FaceAlpha',.3); 
annotation('textarrow',[0.28 0.21],[0.5 0.4])


% PLOTS FOR ALL

% position_cond1 = 1:3:22;
% position_cond2 = 2:3:23;
% position_cond3 = 3:3:24;
% for s=s
%    position_cond1_1 = position_cond1(s);
%    position_cond2_2 = position_cond2(s);
%    position_cond3_3 = position_cond3(s);
% end 
   
% figure(96)
% 
% if c == 1
% subplot(8,3,position_cond1_1)
% if position_cond1_1==1
%    title('Level')
% end
% yyaxis left
% plot(1:length(stride_speed_all(1:420)),stride_speed_all(1:420));
% ylabel([subjs{s}])
% yyaxis right 
% scatter(1:length(STV_OP1.level.STDvalues(:,s)),STV_OP1.level.STDvalues(:,s),'.');
% hold on
% yline(STV_OP1.level.medianvalue(:,s),'r');
% xline(STV_OP1.level.steadystatevalue(:,s));
% % title('COV')
% xlim([0 420])
% end  
% if c == 2
% % subplot(4,1,1)
% subplot(8,3,position_cond2_2)
% if position_cond2_2==2
%    title('Incline')
% end
% yyaxis left
% plot(1:length(stride_speed_all(1:420)),stride_speed_all(1:420));
% yyaxis right 
% scatter(1:length(STV_OP1.incline.STDvalues(:,s)),STV_OP1.incline.STDvalues(:,s),'.');
% hold on
% yline(STV_OP1.incline.medianvalue(:,s),'r');
% xline(STV_OP1.incline.steadystatevalue(:,s));
% xlim([0 420])
% end
% if c == 3 
% subplot(8,3,position_cond3_3)
% if position_cond3_3==3
%    title('Decline')
% end
% yyaxis left
% plot(1:length(stride_speed_all(1:420)),stride_speed_all(1:420));
% yyaxis right 
% scatter(1:length(STV_OP1.decline.STDvalues(:,s)),STV_OP1.decline.STDvalues(:,s),'.');
% hold on
% yline(STV_OP1.decline.medianvalue(:,s),'r');
% xline(STV_OP1.decline.steadystatevalue(:,s));
% xlim([0 420])
% end

%% stride speed plot
if c == 1 
stride_speed_all_barplot_level(s) = mean(stride_speed_all(STV_OP1.level.steadystatevalue(s):(STV_OP1.level.steadystatevalue(s)+400)));
end 
if c == 2 
stride_speed_all_barplot_incline(s) = mean(stride_speed_all(STV_OP1.incline.steadystatevalue(s):(STV_OP1.incline.steadystatevalue(s)+400)));
end 
if c == 3 
stride_speed_all_barplot_decline(s) = mean(stride_speed_all(STV_OP1.decline.steadystatevalue(s):(STV_OP1.decline.steadystatevalue(s)+400)));
end 
%% Stride speeds, lengths, frequency for all subjects 


     if c == 1
     speed_strides.level(s,:) = stride_speed_all(STV_OP1.level.steadystatevalue(s):(STV_OP1.level.steadystatevalue(s)+400));
     end
     if c == 2
     speed_strides.incline(s,:) = stride_speed_all(STV_OP1.incline.steadystatevalue(s):(STV_OP1.incline.steadystatevalue(s)+400));
     end
     if c == 3
     speed_strides.decline(s,:) = stride_speed_all(STV_OP1.decline.steadystatevalue(s):(STV_OP1.decline.steadystatevalue(s)+400));
     end
     
     if c == 1
     length_strides.level(s,:) = stridelength_all(STV_OP1.level.steadystatevalue(s):(STV_OP1.level.steadystatevalue(s)+400));
     end
     if c == 2
     length_strides.incline(s,:) = stridelength_all(STV_OP1.incline.steadystatevalue(s):(STV_OP1.incline.steadystatevalue(s)+400));
     end
     if c == 3
     length_strides.decline(s,:) = stridelength_all(STV_OP1.decline.steadystatevalue(s):(STV_OP1.decline.steadystatevalue(s)+400));
     end
     
     if c == 1
     width_strides.level(s,:) = stridewidth_all(STV_OP1.level.steadystatevalue(s):(STV_OP1.level.steadystatevalue(s)+400));
     end
     if c == 2
     width_strides.incline(s,:) = stridewidth_all(STV_OP1.incline.steadystatevalue(s):(STV_OP1.incline.steadystatevalue(s)+400));
     end
     if c == 3
     width_strides.decline(s,:) = stridewidth_all(STV_OP1.decline.steadystatevalue(s):(STV_OP1.decline.steadystatevalue(s)+400));
     end
     
     
     if c == 1
     time_strides.level(s,:) = stridelength_all(STV_OP1.level.steadystatevalue(s):(STV_OP1.level.steadystatevalue(s)+400));
     end
     if c == 2
     time_strides.incline(s,:) = stridelength_all(STV_OP1.incline.steadystatevalue(s):(STV_OP1.incline.steadystatevalue(s)+400));
     end
     if c == 3
     time_strides.decline(s,:) = stridelength_all(STV_OP1.decline.steadystatevalue(s):(STV_OP1.decline.steadystatevalue(s)+400));
     end
%% BEST FIT MODEL 
% v = 0.5:0.1:1.8;
%     figure(100)
%     hold on
%     if c == 1 
%     subplot(4,2,s)
%     plot(speed_strides.level(s,:),length_strides.level(s,:),'.','DisplayName',[subjs{s} ' ' conds{c}],'Color',color(c,:))
%     hold on
%     [P,fh] = fitsinglemodelprocess_sl(length_strides.level(s,:),speed_strides.level(s,:));
%     plot(v,fh(v,P),'Color',color(c,:));
%     if s == 8 
%     xlabel('Speed')
%     ylabel('Stride length')
%     end
%     title([subjs{s}])
%     end 
%     if c == 2 
%     subplot(4,2,s)
%     plot(speed_strides.incline(s,:),length_strides.incline(s,:),'.','DisplayName',[subjs{s} ' ' conds{c}],'Color',color(c,:))
%     hold on
%     [P,fh] = fitsinglemodelprocess_sl(length_strides.incline(s,:),speed_strides.incline(s,:));
%     plot(v,fh(v,P),'Color',color(c,:));
%     end 
%     if c == 3 
%     subplot(4,2,s)
%     plot(speed_strides.decline(s,:),length_strides.decline(s,:),'.','DisplayName',[subjs{s} ' ' conds{c}],'Color',color(c,:))
%     hold on
%     legend()
%     [P,fh] = fitsinglemodelprocess_sl(length_strides.decline(s,:),speed_strides.decline(s,:));
%     plot(v,fh(v,P),'Color',color(c,:));
%     end 
%     if c == 1 
%     p.level(s,:) = P;
%     end
%      if c == 2 
%     p.incline(s,:) = P;
%      end
%      if c == 3
%     p.decline(s,:) = P;
%      end
%% FINDING VARIABILITY FOR STRIDE LENGTH 
     if c == 1
     [P,fh] = fitsinglemodelprocess_sl(length_strides.level(s,:),speed_strides.level(s,:));
     end 
     if c == 2
     [P,fh] = fitsinglemodelprocess_sl(length_strides.incline(s,:),speed_strides.incline(s,:));
     end 
     if c == 3
     [P,fh] = fitsinglemodelprocess_sl(length_strides.decline(s,:),speed_strides.decline(s,:));
     end 
%     plot(v,fh(v,P),'Color',color(c,:));
    if c == 1 
    p.level(s,:) = P;
    end
     if c == 2 
    p.incline(s,:) = P;
     end
     if c == 3
    p.decline(s,:) = P;
     end
     
     if c == 1
     variation_strides.level.slminusfit(s,:) = var(length_strides.level(s,:) - fh(speed_strides.level(s,:),P));
     variation_strides.level.speedtrend(s,:) = var(fh(speed_strides.level(s,:),P));
     variation_strides.level.totalvar(s,:) = var(length_strides.level(s,:));
     end
     if c == 2
     variation_strides.incline.slminusfit(s,:) = var(length_strides.incline(s,:) - fh(speed_strides.incline(s,:),P));
     variation_strides.incline.speedtrend(s,:) = var(fh(speed_strides.incline(s,:),P));
     variation_strides.incline.totalvar(s,:) = var(length_strides.incline(s,:));
     end
     if c == 3
     variation_strides.decline.slminusfit(s,:) =  var(length_strides.decline(s,:) - fh(speed_strides.decline(s,:),P));
     variation_strides.decline.speedtrend(s,:) = var(fh(speed_strides.decline(s,:),P));
     variation_strides.decline.totalvar(s,:) = var(length_strides.decline(s,:));
     end
%% FINDING VARIABILITY FOR STRIDE WIDTH
     if c == 1
     [W,fh_width] = fitsinglemodelprocess_sw(width_strides.level(s,:),speed_strides.level(s,:));
     end 
     if c == 2
     [W,fh_width] = fitsinglemodelprocess_sw(width_strides.incline(s,:),speed_strides.incline(s,:));
     end 
     if c == 3
     [W,fh_width] = fitsinglemodelprocess_sw(width_strides.decline(s,:),speed_strides.decline(s,:));
     end 
%     plot(v,fh(v,P),'Color',color(c,:));
    if c == 1 
    p.level_width(s,:) = W;
    end
     if c == 2 
    p.incline_width(s,:) = W;
     end
     if c == 3
    p.decline_width(s,:) = W;
     end
     
     if c == 1
     variation_strides_width.level.slminusfit(s,:) = var(width_strides.level(s,:) - fh_width(speed_strides.level(s,:),W));
     variation_strides_width.level.speedtrend(s,:) = var(fh_width(speed_strides.level(s,:),W));
     variation_strides_width.level.totalvar(s,:) = var(width_strides.level(s,:));
     end
     if c == 2
     variation_strides_width.incline.slminusfit(s,:) = var(length_strides.incline(s,:) - fh_width(speed_strides.incline(s,:),W));
     variation_strides_width.incline.speedtrend(s,:) = var(fh_width(speed_strides.incline(s,:),W));
     variation_strides_width.incline.totalvar(s,:) = var(length_strides.incline(s,:));
     end
     if c == 3
     variation_strides_width.decline.slminusfit(s,:) =  var(length_strides.decline(s,:) - fh_width(speed_strides.decline(s,:),W));
     variation_strides_width.decline.speedtrend(s,:) = var(fh_width(speed_strides.decline(s,:),W));
     variation_strides_width.decline.totalvar(s,:) = var(length_strides.decline(s,:));
     end


    end 
end
% 
% fh = @(x,p) p(1)*(x.^p(2));
% % stride_length = sl(:,1)+ sl(:,2);
% % stride_time = mean(steplength_speed,2);
% level_avg_fit = [mean(p.level(:,1)), mean(p.level(:,2))];
% incline_avg_fit = [mean(p.incline(:,1)), mean(p.incline(:,2))];
% decline_avg_fit = [mean(p.decline(:,1)), mean(p.decline(:,2))];
% v = 0.5:0.1:1.8;
% figure(100)
% % subplot(1,2,2)
% plot(v,fh(v,level_avg_fit),'DisplayName',[conds{1}],'Color',color(1,:));
% hold on
% plot(v,fh(v,incline_avg_fit),'DisplayName',[conds{2}],'Color',color(2,:));
% hold on
% plot(v,fh(v,decline_avg_fit),'DisplayName',[conds{3}],'Color',color(3,:));
% legend()
% xlabel('Speed')
% ylabel('Stride length')
% %  
%     
% title('Speed vs Stride length Average')


% fit_table_all = [fit_table_1; fit_table_2];
% fit_table_all = [fit_table_all; fit_table_3];
% lm_all = fitlm(fit_table_all,'Var5~ (Var1 * Var2) + slope','RobustOpts','on')

%% PLOTTING COV BAR GRAPHS
% for i=1:8
% sl_cov_bar_mean(i,:) = [mean(sl_cov_decline{i}); mean(sl_cov_level{i}); mean(sl_cov_incline{i})];
% sl_cov_bar_std(i,:) =  [std(sl_cov_decline{i}); std(sl_cov_level{i}); std(sl_cov_incline{i})];
% sw_cov_bar_mean(i,:) = [mean(sw_cov_decline{i}); mean(sw_cov_level{i}); mean(sw_cov_incline{i})];
% sw_cov_bar_std(i,:) =  [std(sw_cov_decline{i}); std(sw_cov_level{i}); std(sw_cov_incline{i})];


%% STD method BAR plot 
% figure
% plot([STV_OP1.decline.steadystatevalue;STV_OP1.level.steadystatevalue;STV_OP1.incline.steadystatevalue],'-O','LineWidth',2)
% legend(subjs)
% hold on 
% errorbar(mean([STV_OP1.decline.steadystatevalue;STV_OP1.level.steadystatevalue;STV_OP1.incline.steadystatevalue],2),std([STV_OP1.decline.steadystatevalue;STV_OP1.level.steadystatevalue;STV_OP1.incline.steadystatevalue],0,2),'.','Color','k')
% bar(mean([STV_OP1.decline.steadystatevalue;STV_OP1.level.steadystatevalue;STV_OP1.incline.steadystatevalue],2),'FaceColor','none')
% title('Steady State')
% xlabel('Slopes')
% ylabel('Steady State Stride #')
% cov_bar_labels = {'Decline', 'Level', 'Incline'};
% set(gca, 'XTick', 1:3, 'XTickLabel', cov_bar_labels)
%% SPEED BAR PLOT 
% figure
% plot([stride_speed_all_barplot_decline;stride_speed_all_barplot_level;stride_speed_all_barplot_incline],'-O','LineWidth',2)
% legend(subjs)
% hold on 
% errorbar(mean([stride_speed_all_barplot_decline;stride_speed_all_barplot_level;stride_speed_all_barplot_incline],2),std([stride_speed_all_barplot_decline;stride_speed_all_barplot_level;stride_speed_all_barplot_incline],0,2),'.','Color','k')
% bar(mean([stride_speed_all_barplot_decline;stride_speed_all_barplot_level;stride_speed_all_barplot_incline],2),'FaceColor','none')
% title('Speed')
% xlabel('Slopes')
% ylabel('Avg Speed')
% cov_bar_labels = {'Decline', 'Level', 'Incline'};
% set(gca, 'XTick', 1:3, 'XTickLabel', cov_bar_labels)

%% STD method BAR plot for TOTAL NUMBER  

% figure
% plot([STV_OP1.decline.steadystatevalue_all;STV_OP1.level.steadystatevalue_all;STV_OP1.incline.steadystatevalue_all],'-O','LineWidth',2)
% legend(subjs)
% hold on 
% errorbar(mean([STV_OP1.decline.steadystatevalue_all;STV_OP1.level.steadystatevalue_all;STV_OP1.incline.steadystatevalue_all],2),std([STV_OP1.decline.steadystatevalue_all;STV_OP1.level.steadystatevalue_all;STV_OP1.incline.steadystatevalue_all],0,2),'.','Color','k')
% bar(mean([STV_OP1.decline.steadystatevalue_all;STV_OP1.level.steadystatevalue_all;STV_OP1.incline.steadystatevalue_all],2),'FaceColor','none')
% title('Total Steady State Strides')
% xlabel('Slopes')
% ylabel('Percent of Steady Strides')
% cov_bar_labels = {'Decline', 'Level', 'Incline'};
% set(gca, 'XTick', 1:3, 'XTickLabel', cov_bar_labels)


%% VARIANCE PLOTS

% err = [std(variation_strides.decline.totalvar),std(variation_strides.decline.slminusfit),std(variation_strides.decline.speedtrend);std(variation_strides.level.totalvar),std(variation_strides.level.slminusfit),std(variation_strides.level.speedtrend);std(variation_strides.incline.totalvar),std(variation_strides.incline.slminusfit),std(variation_strides.incline.speedtrend)];
% OOOO = [mean(variation_strides.decline.totalvar),mean(variation_strides.decline.slminusfit),mean(variation_strides.decline.speedtrend);mean(variation_strides.level.totalvar),mean(variation_strides.level.slminusfit),mean(variation_strides.level.speedtrend);mean(variation_strides.incline.totalvar),mean(variation_strides.incline.slminusfit),mean(variation_strides.incline.speedtrend)];
% figure 
% bar(OOOO)
% hold on
% ngroups = size(OOOO, 1);
% nbars = size(OOOO, 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, OOOO(:,i), err(:,i), '.');
% end
% hold off
% ylim([0 0.025])
% xticklabels({'Decline','Level','Incline'})

%% STRIDE WIDTH VARIANCE PLOTS
% err = [std(variation_strides_width.decline.totalvar),std(variation_strides_width.decline.slminusfit),std(variation_strides_width.decline.speedtrend);std(variation_strides_width.level.totalvar),std(variation_strides_width.level.slminusfit),std(variation_strides_width.level.speedtrend);std(variation_strides_width.incline.totalvar),std(variation_strides_width.incline.slminusfit),std(variation_strides_width.incline.speedtrend)];
% OOOO = [mean(variation_strides_width.decline.totalvar),mean(variation_strides_width.decline.slminusfit),mean(variation_strides_width.decline.speedtrend);mean(variation_strides_width.level.totalvar),mean(variation_strides_width.level.slminusfit),mean(variation_strides_width.level.speedtrend);mean(variation_strides_width.incline.totalvar),mean(variation_strides_width.incline.slminusfit),mean(variation_strides_width.incline.speedtrend)];
% figure 
% bar(OOOO)
% hold on
% ngroups = size(OOOO, 1);
% nbars = size(OOOO, 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, OOOO(:,i), err(:,i), '.');
% end
% hold off
% ylim([0 0.025])
% xticklabels({'Decline','Level','Incline'})
%% DETRENDED PLOTS 
% figure 
% bar([mean(variation_strides_width.decline.slminusfit);mean(variation_strides_width.level.slminusfit);mean(variation_strides_width.incline.slminusfit)]);
% hold on 
% errorbar([mean(variation_strides_width.decline.slminusfit);mean(variation_strides_width.level.slminusfit);mean(variation_strides_width.incline.slminusfit)],[std(variation_strides_width.decline.slminusfit);std(variation_strides_width.level.slminusfit);std(variation_strides_width.incline.slminusfit)],'.','Color','k');
% ylim([0 0.025])
% xticklabels({'Decline','Level','Incline'})
% 
% figure 
% bar([mean(variation_strides.decline.slminusfit);mean(variation_strides.level.slminusfit);mean(variation_strides.incline.slminusfit)]);
% hold on 
% errorbar([mean(variation_strides.decline.slminusfit);mean(variation_strides.level.slminusfit);mean(variation_strides.incline.slminusfit)],[std(variation_strides.decline.slminusfit);std(variation_strides.level.slminusfit);std(variation_strides.incline.slminusfit)],'.','Color','k');
% hold off
% ylim([0 0.0025])
% xticklabels({'Decline','Level','Incline'})

% figure 
% plot(stride_speed_all,'.')
% xlim([0 450]) 
% ylim([0 1.4])
% figure 
% plot(stridelength_all,'.')
% xlim([0 450])
% ylim([0 1.6])
% figure 
% plot(stridewidth_all,'.')
% xlim([0 450])
% ylim([0 0.3])
%% SELFPACED HIP VS SPEED 
% figure
% plot(datatreadmill.LeftBeltSpeed,(llmarkers.LASI(:,3)+llmarkers.LPSI(:,3)+llmarkers.RASI(:,3)+llmarkers.RPSI(:,3))/4)