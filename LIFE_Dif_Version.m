% close all
clearvars
clc
% % 
% subjs = {'Ll01' 'Ll02' 'Ll03' 'Ll04' 'Ll05' 'Ll06' 'Ll07' 'Ll10'};
subjs = {'Ll04'};
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
        
        conds_notincline = {'level_050' 'level_075'  'level_100' 'level_125' 'level_selfpaced' 'decline_selfpaced' 'decline_075'};
        conds_incline = {'incline_selfpaced' 'incline_075'};
        
      
        %%
        % ML, horizontal/parallel, vertical/perpendicular
        %         F_parallel = forces_df_c.FP2For(:,2)*cosd(theta(c)) + forces_df_c.FP2For(:,3)*sind(theta(c));
        %         F_perpendicular = forces_df_c.FP2For(:,2)*sind(theta(c)) + forces_df_c.FP2For(:,3)*cosd(theta(c));
        %         Rgrf = [forces_df_c.FP2For(:,1) F_parallel F_perpendicular];
        %
        %         F_parallel = forces_df_c.FP1For(:,2)*cosd(theta(c)) + forces_df_c.FP1For(:,3)*sind(theta(c));
        %         F_perpendicular = forces_df_c.FP1For(:,2)*sind(theta(c)) + forces_df_c.FP1For(:,3)*cosd(theta(c));
        %         Lgrf = [forces_df_c.FP1For(:,1) F_parallel F_perpendicular];
        
        %% CUTTING OFF MARKER DATA
%         markers_df.RHEE = markers_df.RHEE(15000:59000,:);
%         markers_df.LHEE = markers_df.LHEE(15000:59000,:);
%         markers_df.RANK = markers_df.RANK(15000:59000,:);
%         markers_df.LANK = markers_df.LANK(15000:59000,:);
%         markers_df.RTOE = markers_df.RTOE(15000:59000,:);
%         markers_df.LTOE = markers_df.LTOE(15000:59000,:);
%         Time_df = Time_df(15000:59000,:);
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
        if Z==sum(strcmp(conds(c),conds_notincline)) 
        [RHS,LTO,LHS,RTO,GE,GEInMiddleDeleted] = GaitEvents_allslopes(Time,llmarkers,markers4GE, HSrefinePre, HSrefinePost, TOminpeakdistance, TOminpeakheight,BW); 
        end 
        
        if Z==sum(strcmp(conds(c),conds_incline)) 
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
      
% h(1)=figure;
% plot(markers_df.RHEE(:,2))
% title('RHEE')
% 
% h(2)=figure;
% plot(markers_df.RTOE(:,2))
% title('RTOE')
% 
% h(3)=figure;
% plot(markers_df.LTOE(:,2))
% title('LTOE')
% 
% h(4)=figure;
% plot(markers_df.LHEE(:,2))
% title('LHEE')
% % 
% savefig(h,'Markers_level_050.fig')            
%% Calc Speed from treadmill

% datatreadmill.Time = interpft(datatreadmill.Time(:,:),length(Frame_df));
% datatreadmill.LeftBeltSpeed = interpft(datatreadmill.LeftBeltSpeed(:,:),length(Frame_df));
% datatreadmill.RightBeltSpeed  = interpft(datatreadmill.RightBeltSpeed(:,:),length(Frame_df));
% datatreadmill.Pitch = interpft(datatreadmill.Pitch(:,:),length(Frame_df));
% 
% datatreadmill.Time = datatreadmill.Time(startidx:stopidx);
% datatreadmill.LeftBeltSpeed =  datatreadmill.LeftBeltSpeed(startidx:stopidx);
% datatreadmill.RightBeltSpeed  = datatreadmill.RightBeltSpeed(startidx:stopidx);
% datatreadmill.Pitch = datatreadmill.Pitch(startidx:stopidx);

% Tread_start_trigger = Time_treadmill(startidx:stopidx);
% Tread_stop_trigger = Time_treadmill)*(stopidx/Total)); 
% 
% datatreadmill.Time = datatreadmill.Time(Tread_start_trigger:Tread_stop_trigger);
% datatreadmill.LeftBeltSpeed = datatreadmill.LeftBeltSpeed(Tread_start_trigger:Tread_stop_trigger);
% datatreadmill.RightBeltSpeed  = datatreadmill.RightBeltSpeed (Tread_start_trigger:Tread_stop_trigger);
% datatreadmill.Pitch = datatreadmill.Pitch(Tread_start_trigger:Tread_stop_trigger);
% % datatreadmill.PitchVGait = datatreadmill.PitchVGait(Tread_start_trigger:Tread_stop_trigger);
% 
% hardstart_treadmill = round(length(datatreadmill.Time)*(Start_Time/Dflow_size));
% hardstop_treadmill = round(length(datatreadmill.Time)*(Stop_Time/Dflow_size));
% 
% datatreadmill.Time = datatreadmill.Time(hardstart_treadmill:hardstop_treadmill);
% datatreadmill.LeftBeltSpeed = datatreadmill.LeftBeltSpeed(hardstart_treadmill:hardstop_treadmill);
% datatreadmill.RightBeltSpeed  = datatreadmill.RightBeltSpeed (hardstart_treadmill:hardstop_treadmill);
% datatreadmill.Pitch = datatreadmill.Pitch(hardstart_treadmill:hardstop_treadmill);

% datatreadmill.PitchVGait = datatreadmill.PitchVGait(hardstart_treadmill:hardstop_treadmill);

% datatreadmill.PitchVGait = interpft(datatreadmill.PitchVGait(:,:),length(Frame_df));
% figure(s)
% hold on 
% scatter(1:length(datatreadmill.RightBeltSpeed(:,1)),datatreadmill.RightBeltSpeed(:,1))
% oo = movmean(datatreadmill.RightBeltSpeed(:,1),240);
% plot(1:length(datatreadmill.RightBeltSpeed(:,1)),oo)

blockSize = [240, 1];
meanFilterFunction = @(theBlockStructure) mean2(theBlockStructure.data(:));
blockAveragedDownSignal = blockproc(datatreadmill.RightBeltSpeed(:,1), blockSize, meanFilterFunction);
[rows, columns] = size(blockAveragedDownSignal);
% hold on
% createfigure(1:length(blockAveragedDownSignal),blockAveragedDownSignal)

% if c == 1
% blockAveragedDownSignal(:,1) = blockAveragedDownSignal(:,1);
% end
% if c == 2
% blockAveragedDownSignal(:,2) = blockAveragedDownSignal(:,1);
% end
% if c == 3
% blockAveragedDownSignal(:,3) = blockAveragedDownSignal(:,1);
% end

% if c == 3 
% figure(s)
% lvlspeed = blockAveragedDownSignal(:,1);
% inclinespeed = blockAveragedDownSignal(:,2);
% declinespeed = blockAveragedDownSignal(:,3);

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

% if c == 2
% customfit = fittype('custom_fit(Ll04_speed(1:end-fitpoints(1),1),0.8,0.1,10,50)');
% f = fit( x, y, customfit);
% end

% [fitresult, gof, output] = createFit1(x1, y1);



% customfit{s,c} = fittwomodelprocess(y1,x1);
% tbl = table(x1,y1,'VariableNames',{'x1','y1'});
% lm = fitlm(tbl);
% [P,lm] = linearreg(y1,x1);

% yf = customfit(1)*(1-exp(-x1/customfit(2)));
% ys = customfit(3)*(1-exp(-x1/customfit(4)));
% 
% figure(1)
% hold on
% plot(x1, yf,'k')
% plot(x1, ys,'g')

% if c == 1
% customfit = custom_fit(Ll04_speed(1:end-fitpoints(1),1),1.0,0.2,8,60);
% end 
% 
% if c == 2
% customfit = custom_fit(Ll04_speed(1:end-fitpoints(1),1),0.8,0.1,10,50);
% end
% 
% if c == 3
% customfit = custom_fit(Ll04_speed(1:end-fitpoints(1),1),0.9,0.2,9,50);
% end
 
% plot(x,customfit,'k')

% [xData, yData] = prepareCurveData( [], y );
% 
% % Set up fittype and options.
% ft = fittype( 'exp2' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% 
% h = plot(fitresult,'k'); %,'DisplayName',[subjs{s} ' ' conds{c} 'fitted'],'Color','k')
% legend()
% xlabel('Time')
% ylabel('Speed m/s')
% title([subjs{s}])

% if c == 3
% savefig(figure(s), [subjs{s} '.fig'])
% end


% lolo = round(length(full_data_treadmill.RightBeltSpeed(1:Stop_Time,1))/10)
% 
% 
% 
% for i = 1:10
%     xx(i) = mean(y(1+lolo*(i-1):lolo*(i)));
%     yy = lolo(i);
% end    
% 
% 
% g = fittype('a-b*exp(-c*x)');
% f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x))]\y; 1]);

% N = length(full_data_treadmill.RightBeltSpeed(1:Stop_Time,1));
% % x = 1:N;
% y = full_data_treadmill.RightBeltSpeed(1:Stop_Time,1);
% p = cell(1, N-60);
% for ix = 1:N-60
%   p{ix} = polyfit(x((0:59)+ix), y((0:59)+ix), 1)';
% end
% p = cell2mat(p)';


% figure(s)
% hold on
% scatter(1:length(full_data_treadmill.RightBeltSpeed(:,1)),full_data_treadmill.RightBeltSpeed(:,1),'DisplayName',[subjs{s} ' ' conds{c}],'MarkerEdgeColor',color(c,:))
% ft=fittype('exp1');
% myfit=fit(x',y,ft)
% plot(myfit,x,y) %color(c,:) (,'DisplayName',[subjs{s} ' ' conds{c} 'fitted'],'Color','k'))
% 
% legend()
% title([subjs{s}])




% f = @(b,x) b(1).*exp(b(2).*x)+b(3);                                     % Objective Function
% B = fminsearch(@(b) norm(y - f(b,x)), [-200; -1; 100])                  % Estimate Parameters
% figure(s)
% hold on
% scatter(1:length(full_data_treadmill.RightBeltSpeed(:,1)),full_data_treadmill.RightBeltSpeed(:,1),'DisplayName',[subjs{s} ' ' conds{c}],'MarkerEdgeColor',color(c,:))
% plot(x, f(B,x), '-r')
% grid
% xlabel('x')
% ylabel('f(x)')
% text(27, 105, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', B))
% title([subjs{s}])

% f = @(b,x) b(1).*exp(b(2).*x)+b(3);                                     % Objective Function
% B = fminsearch(@(b) norm(y - f(b,x)), [-200; -1; 100])

% Estimate Parameters
% figure
% ft=fittype('exp1');
% cf=fit(time,data,ft)
% clear blockAveragedDownSignal
% hold on 
% scatter(1:length(blockAveragedDownSignal(:,1)),inclinespeed)
% scatter(1:length(blockAveragedDownSignal(:,1)),declinespeed)
% 
% title([subjs{s}])
% end
% Tread_start = round(length(datatreadmill.Time)*(Start_Time/Dflow_size));
% Tread_stop = round(length(datatreadmill.Time)*(Stop_Time/Dflow_size));
% 
% datatreadmill.Time = datatreadmill.Time(Tread_start:Tread_stop);
% datatreadmill.LeftBeltSpeed = datatreadmill.LeftBeltSpeed(Tread_start:Tread_stop);
% datatreadmill.RightBeltSpeed  = datatreadmill.RightBeltSpeed (Tread_start:Tread_stop);
% datatreadmill.Pitch = datatreadmill.Pitch(Tread_start:Tread_stop);
% 
% datatreadmill.Time = interpft(datatreadmill.Time(:,:),length(Frame_df));
% datatreadmill.LeftBeltSpeed = interpft(datatreadmill.LeftBeltSpeed(:,:),length(Frame_df));
% datatreadmill.RightBeltSpeed  = interpft(datatreadmill.RightBeltSpeed(:,:),length(Frame_df));
% datatreadmill.Pitch = interpft(datatreadmill.Pitch(:,:),length(Frame_df));
% datatreadmill.PitchVGait = datatreadmill.PitchVGait(Tread_start:Tread_stop);


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
%         if c == 1 || c == 2
%         Time_treadmill(find(datatreadmill.PitchVGait > 0.99*theta(c),1','first'))
%         else
%             Time_treadmill(find(datatreadmill.PitchVGait < 0.99*theta(c),1','first'))
%         end
%         Time_df(end)
%         
%         GEtreadmill = get_GEidx(Time_treadmill, Time_df(GEgood));
%         
%         varnames = fieldnames(datatreadmill);
%         findR = strfind(varnames,'R');
%         for v = 1:length(varnames)
%             if findR{v} == 1, tempidx1 = v; end
%         end
%         
%         findSpeed = strfind(varnames,'Speed');
%         if ~isempty(findSpeed{tempidx1}), rbelt = varnames{tempidx1}; end
%         
%         for sn = 1:length(GEtreadmill)
%             eval(['speeddata = datatreadmill.' rbelt '(GEtreadmill(sn,1):GEtreadmill(sn,5));'])
%             avespeed(sn) = mean(speeddata);
%         end
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
       close(gcf) 
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

%% Significance 
        
        
%         st_halfs(:,1)=st(1:st_length(s,c));
%         st_halfs(:,2)=st(st_length(s,c)+1:end);
%         st_anova(s,c) = anova1(st_halfs);
%         
%         clear st_halfs
%         
%         sp_halfs(:,1)=sp(1:sp_length(s,c));
%         sp_halfs(:,2)=sp(sp_length(s,c)+1:end);
%         sp_anova(s,c) = anova1(sp_halfs);
%         
%         clear sp_halfs
%         
%         sl_halfs(:,1)= sl(1:sl_length(s,c),1);
%         sl_halfs(:,2)= sl(sl_length(s,c)+1:end,1);
%         sl_anova(s,c) = anova1(sl_halfs);
%         
%         clear sl_halfs
%% save figure 

% if c == 1
% h(1)=figure;
% plot(sp(:,1))
% title('Ll02 Level Speed')
% end 
% 
% if c == 2 
% h(2)=figure;
% plot(sp(:,1))
% title('Ll02 Incline Speed')
% end 
% 
% if c == 3 
% h(3)=figure;
% plot(sp(:,1))
% title('Ll02 Decline Speed')
% end 
% 
% 
% savefig(h,'Ll02.fig') 
       
% figure ()
% plot(sp(:,1))
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


% h(s)=figure(s);
% hold on
% plot(sl_all,'.','DisplayName',[subjs{s} ' ' conds{c}],'MarkerEdgeColor',color(c,:));
% xlabel('Steps');
% ylabel('Step Length');
% legend()
% [fitresult, gof, output] = createFit1((1:length(sl_all)),sl_all);
% plot(fitresult);
% title([subjs{s} 'steplength'])

% h(s)=figure(s);
% hold on
% plot(sw_all,'.','DisplayName',[subjs{s} ' ' conds{c}],'MarkerEdgeColor',color(c,:));
% xlabel('Steps');
% ylabel('Step Width');
% legend()
% title([subjs{s} 'Stepwidth'])

% h(s)=figure(s);
% hold on
% plot(st,'.','DisplayName',[subjs{s} ' ' conds{c}],'MarkerEdgeColor',color(c,:));
% xlabel('Steps');
% ylabel('Stride Time');
% legend()
% title([subjs{s} 'stridetime'])

% h(s)=figure(s);
% hold on
% plot(step_time,'.','DisplayName',[subjs{s} ' ' conds{c}],'MarkerEdgeColor',color(c,:));
% xlabel('Steps');
% ylabel('Step Time');
% legend()
% title([subjs{s} 'steptime'])

% h(s)=figure(s);
% hold on
% plot(steplength_speed_all,'.','DisplayName',[subjs{s} ' ' conds{c}],'MarkerEdgeColor',color(c,:));
% xlabel('Steps');
% ylabel('Step Speed');
% legend()
% title([subjs{s} 'StepSpeed'])

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
     

%% PLOTS SPEED VS SL,SW,ST

% Step Speed VS Step Length

% if c == 1;
%     h(1) = figure(1);
%     hold on
%     plot(steplength_speed_all,sl_all,'.','DisplayName',[subjs{s} ' ' conds{c}])
%     legend()
%     [P,fh] = fitsinglemodelprocess(sl_all,steplength_speed_all);
%     plot(steplength_speed_all,fh(steplength_speed_all,P),'DisplayName',[subjs{s} ' ' 'fit']);
%     p{s,c} = sprintf('%0.4f %0.4f',P(1),P(2));
%     title('LEVEL Step Speed vs Step Length')
% end
% fh=[];
% if c == 2;
%     h(2) = figure(2);
%     hold on
%     plot(steplength_speed_all,sl_all,'.','DisplayName',[subjs{s} ' ' conds{c}]);
%     legend()
%     [P,fh] = fitsinglemodelprocess(sl_all,steplength_speed_all);
%     plot(steplength_speed_all,fh(steplength_speed_all,P),'DisplayName',[subjs{s} ' ' 'fit']);
%     p{s,c} = sprintf('%0.4f %0.4f',P(1),P(2));
%     title('INCLINE Step Speed vs Step Length')
% end
% fh=[];
% if c == 3;
%     h(3) = figure(3);
%     hold on
%     plot(steplength_speed_all,sl_all,'.','DisplayName',[subjs{s} ' ' conds{c}]);
%     legend()
%     [P,fh] = fitsinglemodelprocess(sl_all,steplength_speed_all);
%     plot(steplength_speed_all,fh(steplength_speed_all,P),'DisplayName',[subjs{s} ' ' 'fit']);
%     p{s,c} = sprintf('%0.4f %0.4f',P(1),P(2));
%     title('DECLINE Step Speed vs Step Length')
% end
% fh=[];

% Step Speed VS Step Width

% if c == 1;
%     h(1) = figure(1);
%     hold on
%     plot(steplength_speed_all,sw_all,'.','DisplayName',[subjs{s} ' ' conds{c}])
%     legend()
%     [P,fh] = fitsinglemodelprocess(sw_all,steplength_speed_all);
%     plot(steplength_speed_all,fh(steplength_speed_all,P),'DisplayName',[subjs{s} ' ' 'fit']);
%     p{s,c} = sprintf('%0.4f %0.4f',P(1),P(2));
%     title('LEVEL Step Speed vs Step Width')
% end
% if c == 2;
%     h(2) = figure(2);
%     hold on
%     plot(steplength_speed_all,sw_all,'.','DisplayName',[subjs{s} ' ' conds{c}])
%     legend()
%     [P,fh] = fitsinglemodelprocess(           % Estimate Parameters
% figure(s)
% hold on
% scatter(1:length(full_data_treadmill.RightBeltSpeed(:,1)),full_data_treadmill.RightBeltSpeed(:,1),'DisplayName',[subjs{s} ' ' conds{c}],'MarkerEdgeColor',color(c,:))
% plot(x, f(B,x), '-r')
% grid
% xlabel('x')
% ylabel('f(x)')
% text(27, 105, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', B))
% title([subjs{s}])

% f = @(b,x) b(1).*exp(b(2).*x)+b(3);                                     % Objective Function
% B = fminsearch(@(b) norm(y - f(b,x)), [-200; -1; 100])

% Estimate Parameters
% figure
% ft=fittype('exp1');
% cf=fit(time,data,ft)
% clear blockAveragedDownSignal
% hold on 
% scatter(1:length(blockAveragedDownSignal(:,1)),inclinespeed)
% scatter(1:length(blockAveragedDownSignal(:,1)),declinespeed)
% 
% title([subjs{s}])
% end
% Tread_start = round(length(datatreadmill.Time)*(Start_Time/Dflow_size));
% Tread_stop = round(length(datatreadmill.Time)*(Stop_Time/Dflow_size));
% 
% datatreadmill.Time = datatreadmill.Time(Tread_start:Tread_stop);
% datatreadmill.LeftBeltSpeed = datatreadmill.LeftBeltSpeed(Tread_start:Tread_stop);
% datatreadmill.RightBeltSpeed  = datatreadmill.RightBeltSpeed (Tread_start:Tread_stop);
% datatreadmill.Pitch = datatreadmill.Pitch(Tread_start:Tread_stop);
% 
% datatreadmill.Time = interpft(datatreadmill.Time(:,:),length(Frame_df));
% datatreadmill.LeftBeltSpeed = interpft(datatreadmill.LeftBeltSpeed(:,:),length(Frame_df));
% datatreadmill.RightBeltSpeed  = interpft(datatreadmill.RightBeltSpeed(:,:),length(Frame_df));
% datatreadmill.Pitch = interpft(datatreadmill.Pitch(:,:),length(Frame_df));
% datatreadmill.PitchVGait = datatreadmill.PitchVGait(Tread_start:Tread_stop);

% if c == 1 
%     slope = repmat(pitchval(c),[length(sl_all) 1]);
%     fit_table_1 = table(steplength_time_all', sl_all', sw_all', slope, steplength_speed_all');
% %     lm_1 = fitlm(fit_table_1);
% end
% slope = [];
% if c == 2 
%     slope = repmat(pitchval(c),[length(sl_all) 1]);
%     fit_table_2 = table(steplength_time_all', sl_all', sw_all', slope, steplength_speed_all');
% %     lm_2 = fitlm(fit_table_2);
% end
% slope = [];
% if c == 3 
%     slope = repmat(pitchval(c),[length(sl_all) 1]);
%     fit_table_3 = table(steplength_time_all', sl_all', sw_all', slope, steplength_speed_all');
% %     lm_3 = fitlm(fit_table_3);
% end
%% SAVING FILES 

% filename = [subjs(s) conds(c)];
% filename = strjoin(filename,{'_'});
% filename =  strcat(filename,'.mat');
% save(filename,'GEgood')
% GEgood = [];
%% COV

% Treadmill Speed 
if c == 1
  treadmill_speed_cov_level{s} = movcv(full_treadmill_speed.level(:));
end
if c == 2      
  treadmill_speed_cov_incline{s} = movcv(full_treadmill_speed.incline(:));
end
if c == 3
  treadmill_speed_cov_decline{s} = movcv(full_treadmill_speed.decline(:));
end


% Step Speed 
% window_cov_speed= floor(length(steplength_speed_all)/50);
% for i=(1:50)
%     if c == 1
%         if i==1
%             speed_cov_level(i,s) = getCV(steplength_speed_all(1:window_cov_speed));
%         else
%             speed_cov_level(i,s) = getCV(steplength_speed_all(((i-1)*window_cov_speed):(i*window_cov_speed)));
%         end
%     end
%     if c == 2
%         if i==1
%             speed_cov_incline(i,s) = getCV(steplength_speed_all(1:window_cov_speed));
%         else
%             speed_cov_incline(i,s) = getCV(steplength_speed_all(((i-1)*window_cov_speed):(i*window_cov_speed)));
%         end
%     end
%     if c == 3
%         if i==1
%             speed_cov_decline(i,s) = getCV(steplength_speed_all(1:window_cov_speed));
%         else
%             speed_cov_decline(i,s) = getCV(steplength_speed_all(((i-1)*window_cov_speed):(i*window_cov_speed)));
%         end
%     end
% end 

% Step Speed 
if c == 1
  speed_cov_level{s} = movcv(steplength_speed_all(:));
end
if c == 2      
  speed_cov_incline{s} = movcv(steplength_speed_all(:));      
end
if c == 3
  speed_cov_decline{s} = movcv(steplength_speed_all(:));        
end


% Step Length 
if c == 1
  sl_cov_level{s} = movcv(sl_all(:));
end
if c == 2      
  sl_cov_incline{s} = movcv(sl_all(:));      
end
if c == 3
  sl_cov_decline{s} = movcv(sl_all(:));        
end


% Step Width 
if c == 1
  sw_cov_level{s} = movcv(sw_all(:));
end
if c == 2      
  sw_cov_incline{s} = movcv(sw_all(:));      
end
if c == 3
  sw_cov_decline{s} = movcv(sw_all(:));        
end


% Step Time
if c == 1
  steptime_cov_level{s} = movcv(steplength_time_all(:));
end
if c == 2      
  steptime_cov_incline{s} = movcv(steplength_time_all(:));      
end
if c == 3
  steptime_cov_decline{s} = movcv(steplength_time_all(:));        
end



%  
% figure(100)
% sgtitle('Coeff of Var')
% 
% % Plot for treadmill speed 
% % if c == 1
% % subplot(5,3,1)
% % hold on
% % scatter(1:length(treadmill_speed_cov_level{s}),treadmill_speed_cov_level{s},'filled')
% % ylabel('Treadmill Speed')
% % title('Level')
% % end  
% % if c == 2
% % subplot(5,3,2)
% % hold on
% % scatter(1:length(treadmill_speed_cov_incline{s}),treadmill_speed_cov_incline{s},'filled')
% % title('Incline')
% % end
% % if c == 3 
% % subplot(5,3,3)
% % hold on
% % scatter(1:length(treadmill_speed_cov_decline{s}),treadmill_speed_cov_decline{s},'filled')
% % title('Decline')
% % end
% 
% % Plot for Step Speed 
% if c == 1
% subplot(4,3,1)
% xlim([0 440])
% hold on
% scatter(1:length(speed_cov_level{s}),speed_cov_level{s},'filled')
% ylabel('Step Speed')
% title('Level')
% end  
% if c == 2
% subplot(4,3,2)
% xlim([0 440])
% hold on
% scatter(1:length(speed_cov_incline{s}),speed_cov_incline{s},'filled')
% title('Incline')
% end
% if c == 3 
% subplot(4,3,3)
% xlim([0 440])
% hold on
% scatter(1:length(speed_cov_decline{s}),speed_cov_decline{s},'filled')
% title('Incline')
% end
% 
% % Plot for step length  
% if c == 1
% subplot(4,3,4)
% xlim([0 440])
% hold on
% scatter(1:length(sl_cov_level{s}),sl_cov_level{s},'filled')
% ylabel('Step Length')
% end  
% if c == 2
% subplot(4,3,5)
% xlim([0 440])
% hold on
% scatter(1:length(sl_cov_incline{s}),sl_cov_incline{s},'filled')
% end
% if c == 3 
% subplot(4,3,6)
% xlim([0 440])
% hold on
% scatter(1:length(sl_cov_decline{s}),sl_cov_decline{s},'filled')
% end
% 
% % Plot for step Width   
% if c == 1
% subplot(4,3,7)
% xlim([0 440])
% hold on
% scatter(1:length(sw_cov_level{s}),sw_cov_level{s},'filled')
% ylabel('Step Width')
% end  
% if c == 2
% subplot(4,3,8)
% xlim([0 440])
% hold on
% scatter(1:length(sw_cov_incline{s}),sw_cov_incline{s},'filled')
% end
% if c == 3 
% subplot(4,3,9)
% xlim([0 440])
% hold on
% scatter(1:length(sw_cov_decline{s}),sw_cov_decline{s},'filled')
% end
% 
% % Plot for step time     
% if c == 1
% subplot(4,3,10)
% xlim([0 440])
% hold on
% scatter(1:length(steptime_cov_level{s}),steptime_cov_level{s},'filled')
% ylabel('Step Time')
% xlabel('Time')
% end  
% if c == 2
% subplot(4,3,11)
% xlim([0 440])
% hold on
% scatter(1:length(steptime_cov_incline{s}),steptime_cov_incline{s},'filled')
% end
% if c == 3 
% subplot(4,3,12)
% xlim([0 440])
% hold on
% scatter(1:length(steptime_cov_decline{s}),steptime_cov_decline{s},'filled')
% end
%% COV for individual subjects - FOR ONE SUBJECT AT A TIME
% 
% figure(100)
% sgtitle([subjs{s}])
% 
% % % Plot for treadmill speed 
% % if c == 1
% % subplot(5,1,1)
% % hold on
% % scatter(treadmill_speed_cov_level(:,s))
% % ylabel('Treadmill Speed')
% % % title('Level')
% % end  
% % if c == 2
% % subplot(5,1,1)
% % hold on
% % scatter(treadmill_speed_cov_incline(:,s))
% % % title('Incline')
% % end
% % if c == 3 
% % subplot(5,1,1)
% % hold on
% % scatter(treadmill_speed_cov_decline(:,s))
% % % title('Decline')
% % end
% 
% Plot for Step Speed 
% if c == 1
% subplot(4,1,1)
% hold on
% scatter(1:length(speed_cov_level{s}),speed_cov_level{s},'filled')
% ylabel('Step Speed')
% end  
% if c == 2
% subplot(4,1,1)
% hold on
% scatter(1:length(speed_cov_incline{s}),speed_cov_incline{s},'filled')
% end
% if c == 3 
% subplot(4,1,1)
% hold on
% scatter(1:length(speed_cov_decline{s}),speed_cov_decline{s},'filled')
% end

% 
% % Plot for step length  
% if c == 1
% subplot(4,1,2)
% hold on
% scatter(1:length(sl_cov_level{s}),sl_cov_level{s},'filled')
% ylabel('Step Length')
% end  
% if c == 2
% subplot(4,1,2)
% hold on
% scatter(1:length(sl_cov_incline{s}),sl_cov_incline{s},'filled')
% end
% if c == 3 
% subplot(4,1,2)
% hold on
% scatter(1:length(sl_cov_decline{s}),sl_cov_decline{s},'filled')
% end
% 
% % Plot for step Width   
% if c == 1
% subplot(4,1,3)
% hold on
% scatter(1:length(sw_cov_level{s}),sw_cov_level{s},'filled')
% ylabel('Step Width')
% end  
% if c == 2
% subplot(4,1,3)
% hold on
% scatter(1:length(sw_cov_incline{s}),sw_cov_incline{s},'filled')
% end
% if c == 3 
% subplot(4,1,3)
% hold on
% scatter(1:length(sw_cov_decline{s}),sw_cov_decline{s},'filled')
% end
% 
% % Plot for step time     
% if c == 1
% subplot(4,1,4)
% hold on
% scatter(1:length(steptime_cov_level{s}),steptime_cov_level{s},'filled')
% ylabel('Step Time')
% end  
% if c == 2
% subplot(4,1,4)
% hold on
% scatter(1:length(steptime_cov_incline{s}),steptime_cov_incline{s},'filled')
% end
% if c == 3 
% subplot(4,1,4)
% hold on
% scatter(1:length(steptime_cov_decline{s}),steptime_cov_decline{s},'filled')
% end
% legend('level','incline','decline')
% 
% % % close(gcf) 
%% ALL SUBJECTS COV VS REG
% position_sub_cov = [1 3 5 7 9 11 13 15];
% position_sub_st = [2 4 6 8 10 12 14 16];
% for s=s
%    position_sub_cov_1 = position_sub_cov(s);
%    position_sub_st_1 = position_sub_st(s);
% end 
%    
% figure(98)    
% if c == 1
% % subplot(4,1,1)
% subplot(8,2,position_sub_cov_1)
% if position_sub_cov_1==1
%    title('CV (100)')
% end
% hold on
% scatter(1:length(speed_cov_level{s}), speed_cov_level{s},'.')
% ylabel([subjs{s}])
% % title('COV')
% xlim([0 400])
% end  
% if c == 2
% % subplot(4,1,1)
% subplot(8,2,position_sub_cov_1)
% hold on
% scatter(1:length(speed_cov_incline{s}),speed_cov_incline{s},'.')
% xlim([0 400])
% end
% if c == 3 
% subplot(8,2,position_sub_cov_1)
% hold on
% scatter(1:length(speed_cov_decline{s}),speed_cov_decline{s},'.')
% xlim([0 400])
% end
% legend('level','incline','decline')
% 
% if c == 1
% % subplot(4,1,1)
% subplot(8,2,position_sub_st_1)
% if position_sub_cov_1==1
%    title('Actual Step Speed')
% end
% hold on
% scatter(1:length((steplength_speed_all(:))),(steplength_speed_all(:)),'.')
% % ylabel('Step Speed')
% % title('Step Length')
% xlim([0 400])
% end  
% if c == 2
% % subplot(4,1,1)
% subplot(8,2,position_sub_st_1)
% hold on
% scatter(1:length((steplength_speed_all(:))),(steplength_speed_all(:)),'.')
% xlim([0 400])
% end
% if c == 3 
% subplot(8,2,position_sub_st_1)
% hold on
% scatter(1:length((steplength_speed_all(:))),(steplength_speed_all(:)),'.')
% xlim([0 400])
% end
%% HISTOGRAM
% sgtitle('Step Time (0.05)')
% edges = [0.35:0.05:1.65];
% % figure(c)
% figure(97)
% if c == 1
% subplot(8,3, start_location)
% h.level.subjs(s) = histogram(steplength_speed_all,edges);
% ylim([0 210])
% ylabel([subjs{s}])
% end
% % figure(c)
% if c == 2
% subplot(8,3, start_location)
% h.incline.subjs(s) = histogram(steplength_speed_all,edges);
% ylim([0 210])
% end
% % figure(c)
% if c == 3
% subplot(8,3, start_location)
% h.decline.subjs(s) = histogram(steplength_speed_all,edges);
% ylim([0 210])
% end
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

     
%% all plots for first and second half speeds 
% figure(96)
% 
% if c == 1
% subplot(8,3,position_cond1_1)
% plot(speed_strides.level(s,(1:200)),'k');
% hold on 
% plot(speed_strides.level(s,(201:400)),'g');
% if position_cond1_1==1
%    title('Level')
% end
% ylabel([subjs{s}])
% % title('COV')
% xlim([0 200])
% end  
% if c == 2
% % subplot(4,1,1)
% subplot(8,3,position_cond2_2)
% plot(speed_strides.incline(s,(1:200)),'k');
% hold on 
% plot(speed_strides.incline(s,(201:400)),'g');
% if position_cond2_2==2
%    title('Incline')
% end
% ylabel([subjs{s}])
% % title('COV')
% xlim([0 200])
% end
% if c == 3 
% subplot(8,3,position_cond3_3)
% plot(speed_strides.decline(s,(1:200)),'k');
% hold on 
% plot(speed_strides.decline(s,(201:400)),'g');
% if position_cond3_3==3
%    title('Decline')
% end
% ylabel([subjs{s}])
% % title('COV')
% xlim([0 200])
% end
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

%% significance half speeds 


% for i = 1:8
% ttt_1(i) = ttest(time_strides.level(i,(1:200)),time_strides.incline(i,(201:400)));
% end 
% for i = 1:8
% ttt_2(i) = ttest(time_strides.incline(i,(1:200)),time_strides.level(i,(201:400)));
% end 
% for i = 1:8
% ttt_3(i) = ttest(time_strides.decline(i,(1:200)),time_strides.decline(i,(201:400)));
% end 
% 
% for i = 1:8
% ttt_1_1sthalf(i) = ttest(time_strides.level(i,(1:200)));
% ttt_1_2ndhalf(i) = ttest(time_strides.level(i,(201:400)));
% end 
% for i = 1:8
% ttt_2_1sthalf(i) = ttest(time_strides.incline(i,(1:200)));
% ttt_2_2ndhalf(i) = ttest(time_strides.incline(i,(201:400)));
% end 
% for i = 1:8
% ttt_3_1sthalf(i) = ttest(time_strides.decline(i,(1:200)));
% ttt_3_2ndhalf(i) = ttest(time_strides.decline(i,(201:400)));
% end 
% 
% for i = 1:8
% ttt_1_ONE(i) = ttest(time_strides.level(i,:));
% end 
% for i = 1:8
% ttt_2_ONE(i) = ttest(time_strides.incline(i,:));
% end 
% for i = 1:8
% ttt_3_ONE(i) = ttest(time_strides.decline(i,:));
% end 
%% plots half speed 
% for i = 1:8
% figure 
% bar(mean(speed_strides.level(i,(1:200))),'EdgeColor','k','FaceColor','none');
% hold on 
% bar(mean(speed_strides.level(i,(201:400))),'EdgeColor','g','FaceColor','none');
% end 
% for i = 1:8
% figure 
% bar(mean(speed_strides.incline(i,(1:200))),'EdgeColor','k','FaceColor','none');
% hold on 
% bar(mean(speed_strides.incline(i,(201:400))),'EdgeColor','g','FaceColor','none');
% end 
% for i = 1:8
% figure
% bar(mean(speed_strides.decline(i,(1:200))),'EdgeColor','k','FaceColor','none');
% hold on 
% bar(mean(speed_strides.decline(i,(201:400))),'EdgeColor','g','FaceColor','none');
% end

% figure
% f = 1:3:24;
% for i = 1:8
% subplot(8,3,f(i))
% b=bar([mean(speed_strides.level(1,(1:200))) mean(speed_strides.level(1,(201:400)))]);
% b.FaceColor = 'flat';
% b.CData(:,2) = [.5 0];
% ylim([0 1.5])
% xticklabels({'first','second'})
% if i == 1 
%     title('Level')
% end
% end 
% for i = 1:8
% subplot(8,3,f(i)+1) 
% b=bar([mean(speed_strides.incline(i,(1:200))) mean(speed_strides.incline(i,(201:400)))]);
% b.FaceColor = 'flat';
% b.CData(:,2) = [.5 0];
% ylim([0 1.5])
% xticklabels({'first','second'})
% if i == 1 
%     title('Incline')
% end
% end 
% for i = 1:8
% subplot(8,3,f(i)+2)
% b=bar([mean(speed_strides.decline(i,(1:200))) mean(speed_strides.decline(i,(201:400)))]);
% b.FaceColor = 'flat';
% b.CData(:,2) = [.5 0];
% ylim([0 1.5])
% xticklabels({'first','second'})
% if i == 1 
%     title('Decline')
% end
% end

% for i = 1:8
% figure 
% plot(speed_strides.level(i,(1:200)),'k');
% hold on 
% plot(speed_strides.level(i,(201:400)),'g');
% end 
% for i = 1:8
% figure 
% plot(speed_strides.incline(i,(1:200)),'k');
% hold on 
% plot(speed_strides.incline(i,(201:400)),'g');
% end 
% for i = 1:8
% figure
% plot(speed_strides.decline(i,(1:200)),'k');
% hold on 
% plot(speed_strides.decline(i,(201:400)),'g');
% end

% for i = 1:8
% mean(speed_strides.level(i,(1:200)))
% hold on 
% plot(speed_strides.level(i,(201:400)),'g');
% end 
% for i = 1:8
% figure 
% plot(speed_strides.incline(i,(1:200)),'k');
% hold on 
% plot(speed_strides.incline(i,(201:400)),'g');
% end 
% for i = 1:8
% figure
% plot(speed_strides.decline(i,(1:200)),'k');
% hold on 
% plot(speed_strides.decline(i,(201:400)),'g');
% end


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