%% SCRIPT FOR BME6215 to plot walking kinematics & GRF
%
% This script reads in D-flow data and animates a mocap skeleton.
%
% IMPORTANT, in the BRaIN Lab
% if walking towards the front of the treadmill
% x = side to side, mediolateral
% y = vertical
% z = forwards and backwards, anterior-posterior
%
% (+x, right; -x, left)
% (+y up; -y down)
% (+z backwards; -z forwards)

% This script uses code sections, created with the %%
% In the highlighted section, you can execute just that section with
% cmd+enter (on mac)
%
% created by Helen J. Huang, 03/31/17

% close all
clearvars
clc

subjs = {'Ll01' 'Ll02' 'Ll03' 'Ll04' 'Ll06' 'Ll07' 'Ll10'};
conds = {'level_050' 'level_075' 'level_100' 'level_125' 'level_selfpaced' 'incline_selfpaced' 'decline_selfpaced' 'incline_075' 'decline_075'};
theta = [-9 0 9];
projfolder = pwd;

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
        [Frame_df, Time_df, markers_df, forces_df] = import_dflow(dflow_file);
        
        datatreadmill = import_dflow_treadmill(dflow_treadmill_file);
        Time_treadmill = datatreadmill.Time-datatreadmill.Time(1);
        
        %% Convert to conventional coords
        for m = 1:length(markers_df.labels)
            eval(['markers_df_c.' markers_df.labels{m} ' = convert_coords2conventional(markers_df. ' markers_df.labels{m} ');']);
        end
        
        for m = 1:length(forces_df.labels)
            eval(['forces_df_c.' forces_df.labels{m} ' = convert_coords2conventional(forces_df. ' forces_df.labels{m} ');']);
        end
        
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

        [GE.dfgrf]=get_gaitEvents_GRF(forces_df.FP1For,forces_df.FP2For);
        
        RHS = GE.dfgrf(:,1);
        LTO = GE.dfgrf(:,2);
        LHS = GE.dfgrf(:,3);
        RTO = GE.dfgrf(:,4);
        
        GEgood = gaiteventCheck3(RHS, LTO, LHS, RTO, forces_df.FP2For(:,2), forces_df.FP1For(:,2));
       
        acceptable = GEgood(:,1)>7200;
        GEgood = GEgood(acceptable(:),:);
%         if GEgood(:,1) >= 7200
            
            
        
%         [RHS,LHS,RTO,LTO] = findEventIndicesForWalkingOnTreadmill(forces_df_c.FP2For, forces_df_c.FP1For, 240);
%         
%         if ~isempty(RHS)
%             GEgood = gaiteventCheck3(RHS, LTO, LHS, RTO, forces_df_c.FP2For(:,3), forces_df_c.FP1For(:,3));
%         else
%             [GE.dfgrf]=get_gaitEvents_GRF(forces_df.FP1For,forces_df.FP2For);
%             RHS = GE.dfgrf(:,1);
%             LTO = GE.dfgrf(:,2);
%             LHS = GE.dfgrf(:,3);
%             RTO = GE.dfgrf(:,4);
%             GEgood = gaiteventCheck3(RHS, LTO, LHS, RTO, forces_df.FP2For(:,2), forces_df.FP1For(:,2)); 
%         end
        
        %%
%         if c == 1 || c == 2
%         Time_treadmill(find(datatreadmill.PitchVGait > 0.99*theta(c),1','first'))
%         else
%             Time_treadmill(find(datatreadmill.PitchVGait < 0.99*theta(c),1','first'))
%         end
%         Time_df(end)
%         
        GEtreadmill = get_GEidx(Time_treadmill, Time_df(GEgood));
        
        varnames = fieldnames(datatreadmill);
        findR = strfind(varnames,'R');
        for v = 1:length(varnames)
            if findR{v} == 1, tempidx1 = v; end
        end
        
        findSpeed = strfind(varnames,'Speed');
        if ~isempty(findSpeed{tempidx1}), rbelt = varnames{tempidx1}; end
        
        for sn = 1:length(GEtreadmill)
            eval(['speeddata = datatreadmill.' rbelt '(GEtreadmill(sn,1):GEtreadmill(sn,5));'])
            avespeed(sn) = mean(speeddata);
        end
        %%
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
        
%         close(gcf)
        %%
        
        %% Step width & Step length
        
        clear sl sw badidx1 badidx2 allbadsl
        
%         sl(:,1) = steplength(markers_df_c.RHEE(:,2), markers_df_c.LHEE(:,2), GEgood(:,1));
%         sw(:,1) = stepwidth(markers_df_c.RHEE(:,1), markers_df_c.LHEE(:,1), GEgood(:,1));
%         
%         sl(:,2) = steplength(markers_df_c.LHEE(:,2), markers_df_c.RHEE(:,2), GEgood(:,3));
%         sw(:,2) = stepwidth(markers_df_c.LHEE(:,1), markers_df_c.RHEE(:,1), GEgood(:,3));

        sl(:,1) = steplength(markers_df_c.RHEE(:,2), markers_df_c.LHEE(:,2), GEgood(:,1));
        sw(:,1) = stepwidth(markers_df_c.RHEE(:,1), markers_df_c.LHEE(:,1), GEgood(:,1));
        
        sl(:,2) = steplength(markers_df_c.LHEE(:,2), markers_df_c.RHEE(:,2), GEgood(:,3));
        sw(:,2) = stepwidth(markers_df_c.LHEE(:,1), markers_df_c.RHEE(:,1), GEgood(:,3));
        goodidx1 = find(sl(:,1) > 0.9*mean(sl(:,1)));
        goodidx2 = find(sl(:,2) > 0.9*mean(sl(:,2)));
        goodsl = intersect(goodidx1, goodidx2); 
        tempsl = sl(goodsl,:);
        tempsw = sw(goodsl,:);
        
        sl = tempsl;
        sw = tempsw;
        sp = avespeed(goodsl(:))';
%         sp = avespeed(39:end)';
        for sn = 1:length(GEgood)
            st(sn,1) = Time_df(GEgood(sn,5))-Time_df(GEgood(sn,1));
        end
%         for sn = 41:length(GEgood)
%             st(sn,1) = Time_df(GEgood(sn,5))-Time_df(GEgood(sn,1));
%         end
%         
        
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
        %% Avg Speed
        avgspeed1(:,c)= nanmean(GRP.sp(:,c))
        avgspeed= round(avgspeed1,2)
    end
end 
%%
for c = 1:length(conds)
    aveSL(:,c) = nanmean([squeeze(GRP.sl(:,1,c)) squeeze(GRP.sl(:,2,c))],2);
    aveSW(:,c) = nanmean([squeeze(GRP.sw(:,1,c)) squeeze(GRP.sw(:,2,c))],2);
end

color2use = {'r' 'g' 'c' 'm' 'y' 'k' 'k'};
figure
subplot(121)
hold on
bar(nanmean(aveSL), 'b', 'linewidth', 2)
errorbar(nanmean(aveSL),nanstd(aveSL), 'b','linewidth', 2, 'linestyle','none')
for s = 1:length(subjs)
    plot(aveSL(s,:))
end
ylimits = ylim(gca);
axis([0 4 ylimits(1) ylimits(2)])
title('SL')

subplot(122)
hold on
bar(nanmean(aveSW), 'b', 'linewidth', 2)
errorbar(nanmean(aveSW),nanstd(aveSW), 'b','linewidth', 2, 'linestyle','none')
for s = 1:length(subjs)
    plot(aveSW(s,:))
end
ylimits = ylim(gca);
axis([0 4 ylimits(1) ylimits(2)])
title('SW')

%%
for c = 1:length(conds)
    aveSLvar(:,c) = nanmean([squeeze(GRP.slvar(:,1,c)) squeeze(GRP.slvar(:,2,c))],2);
    aveSWvar(:,c) = nanmean([squeeze(GRP.swvar(:,1,c)) squeeze(GRP.swvar(:,2,c))],2);
end

color2use = {'r' 'g' 'c' 'm' 'r:' 'k' 'k:'};
figure
subplot(121)
hold on
bar(nanmean(aveSLvar), 'b')
errorbar(nanmean(aveSLvar),nanstd(aveSLvar), 'b','linewidth', 2, 'linestyle','none')
for s = 1:length(subjs)
    plot(aveSLvar(s,:))
end
ylimits = ylim(gca);
axis([0 4 ylimits(1) ylimits(2)])
title('SL variability')

subplot(122)
hold on
bar(nanmean(aveSWvar), 'b')
errorbar(nanmean(aveSWvar),nanstd(aveSWvar), 'b', 'linewidth', 2, 'linestyle','none')
for s = 1:length(subjs)
    plot(aveSWvar(s,:))
end
ylimits = ylim(gca);
axis([0 4 ylimits(1) ylimits(2)])
title('SW variability')
%% main plots
% 
% figure
% % errorbar(nanmean(aveSL(:,1:7)),nanstd(aveSL(:,1:7)),'Marker','*','linewidth', 2, 'linestyle','none'),hold on
% % xticklabels({'0.50','0.75','0.93','1.00','1.05','1.23','1.25'})
% % xticklabels({'1.23','0.93','1.05'})
% y=nanmean(aveSL(1:7,1:4));
% x=avgspeed(1:4);
% w=nanmean(aveSL(1:7,5:7));
% z=avgspeed(5:7);
% % a=nanmean(aveSL(:,8:9));
% % b=avgspeed(8:9);
% errorbar(x,y,nanstd(aveSL(1:7,1:4)),'Marker','*','linewidth', 2, 'linestyle','none'),hold on
% errorbar(z,w,nanstd(aveSL(1:7,5:7)),'Marker','*','linewidth', 2, 'linestyle','none'), hold on
% % scatter(avgspeed(1:4),nanmean(aveSL(8:9,1:4)))
% % errorbar(b,a,nanstd(aveSL(:,8:9)),'Marker','*','linewidth', 2, 'linestyle','none')
% xlim([0.4 1.3])
% legend('level-fixedspeed','slope-selfpace','slope-fixedspeed')
% xlabel('speed(m/s)')
% ylabel('Step Length(m)')
% title('SL')
% 
% figure
% % errorbar(nanmean(aveSLvar(:,1:4)),nanstd(aveSLvar(:,1:4)),'Marker','*','linewidth', 2, 'linestyle','none')
% % ylim([0.005 0.03])
% % xlim([0 5])
% % xticks([1 2 3 4])
% % xticklabels({'0.50','0.75','1.00','1.25'})
% y=nanmean(aveSLvar(:,1:4));
% x=avgspeed(1:4);
% w=nanmean(aveSLvar(:,5:7));
% z=avgspeed(5:7);
% % a=nanmean(aveSLvar(:,8:9));
% % b=avgspeed(8:9);
% errorbar(x,y,nanstd(aveSLvar(:,1:4)),'Marker','*','linewidth', 2, 'linestyle','none'),hold on
% errorbar(z,w,nanstd(aveSLvar(:,5:7)),'Marker','*','linewidth', 2, 'linestyle','none'), hold on
% % errorbar(b,a,nanstd(aveSLvar(:,8:9)),'Marker','*','linewidth', 2, 'linestyle','none')
% xlim([0.4 1.3])
% legend('level-fixedspeed','slope-selfpace','slope-fixedspeed')
% xlabel('speed(m/s)')
% ylabel('Step Length variability(m)')
% title('SL Var')
% 
% figure
% % errorbar(nanmean(aveSW(:,1:4)),nanstd(aveSW(:,1:4)),'Marker','*','linewidth', 2, 'linestyle','none')
% % ylim([0 0.3])
% % xlim([0 5])
% % xticks([1 2 3 4])
% % xticklabels({'0.50','0.75','1.00','1.25'})
% y=nanmean(aveSW(:,1:4));
% x=avgspeed(1:4);
% w=nanmean(aveSW(:,5:7));
% z=avgspeed(5:7);
% % a=nanmean(aveSW(:,8:9));
% % b=avgspeed(8:9);
% errorbar(x,y,nanstd(aveSW(:,1:4)),'Marker','*','linewidth', 2, 'linestyle','none'),hold on
% errorbar(z,w,nanstd(aveSW(:,5:7)),'Marker','*','linewidth', 2, 'linestyle','none'), hold on
% % errorbar(b,a,nanstd(aveSW(:,8:9)),'Marker','*','linewidth', 2, 'linestyle','none')
% ylim([0 0.3])
% xlim([0.4 1.3])
% legend('level-fixedspeed','slope-selfpace','slope-fixedspeed')
% xlabel('speed(m/s)')
% ylabel('Step Width(m)')
% title('SW')
% 
% figure
% % errorbar(nanmean(aveSWvar(:,1:4)),nanstd(aveSWvar(:,1:4)),'Marker','*','linewidth', 2, 'linestyle','none')
% % ylim([0 0.04])
% % xlim([0 5])
% % xticks([1 2 3 4])
% % xticklabels({'0.50','0.75','1.00','1.25'})
% y=nanmean(aveSWvar(:,1:4));
% x=avgspeed(1:4);
% w=nanmean(aveSWvar(:,5:7));
% z=avgspeed(5:7);
% % a=nanmean(aveSWvar(:,8:9));
% % b=avgspeed(8:9);
% errorbar(x,y,nanstd(aveSWvar(:,1:4)),'Marker','*','linewidth', 2, 'linestyle','none'),hold on
% errorbar(z,w,nanstd(aveSWvar(:,5:7)),'Marker','*','linewidth', 2, 'linestyle','none'), hold on
% % errorbar(b,a,nanstd(aveSWvar(:,8:9)),'Marker','*','linewidth', 2, 'linestyle','none')
% ylim([0 0.04])
% xlim([0.4 1.3])
% legend('level-fixedspeed','slope-selfpace','slope-fixedspeed')
% xlabel('speed(m/s)')
% ylabel('Step Width Variability(m)')
% title('SW Var')
% 
% figure
% y=nanmean(GRP.st(:,1:4));
% x=avgspeed(1:4);
% w=nanmean(GRP.st(:,5:7));
% z=avgspeed(5:7);
% % a=nanmean(GRP.st(:,8:9));
% % b=avgspeed(8:9);
% errorbar(x,y,nanstd(GRP.st(:,1:4)),'Marker','*','linewidth', 2, 'linestyle','none'),hold on
% errorbar(z,w,nanstd(GRP.st(:,5:7)),'Marker','*','linewidth', 2, 'linestyle','none'), hold on
% % errorbar(b,a,nanstd(GRP.st(:,8:9)),'Marker','*','linewidth', 2, 'linestyle','none')
% ylim([0.9 1.7])
% xlim([0.4 1.3])
% legend('level-fixedspeed','slope-selfpace','slope-fixedspeed')
% xlabel('speed(m/s)')
% ylabel('Stride Period(s)')
% title('SP')
% 
% figure
% y=nanmean(GRP.stvar(:,1:4));
% x=avgspeed(1:4);
% w=nanmean(GRP.stvar(:,5:7));
% z=avgspeed(5:7);
% % a=nanmean(GRP.stvar(:,8:9));
% % b=avgspeed(8:9);
% errorbar(x,y,nanstd(GRP.stvar(:,1:4)),'Marker','*','linewidth', 2, 'linestyle','none'),hold on
% errorbar(z,w,nanstd(GRP.stvar(:,5:7)),'Marker','*','linewidth', 2, 'linestyle','none'), hold on
% % errorbar(b,a,nanstd(GRP.stvar(:,8:9)),'Marker','*','linewidth', 2, 'linestyle','none')
% ylim([0 0.4])
% xlim([0.4 1.3])
% legend('level-fixedspeed','slope-selfpace','slope-fixedspeed')
% xlabel('speed(m/s)')
% ylabel('Stride Period Variability(s)')
% title('SP Variability')
% %% bar graphs
% figure
% SLbaryoung = [nanmean(aveSL(:,2)),nanmean(aveSL(:,5));nanmean(aveSL(:,8)),nanmean(aveSL(:,6));nanmean(aveSL(:,9)),nanmean(aveSL(:,7))];
% bar(SLbaryoung)
% xticklabels({'Level','Incline','Decline'})
% legend('0.75 m/s','selfpace')
% title('SL')
% 
% figure
% SLvarbaryoung = [nanmean(aveSLvar(:,2)),nanmean(aveSLvar(:,5));nanmean(aveSLvar(:,8)),nanmean(aveSLvar(:,6));nanmean(aveSLvar(:,9)),nanmean(aveSLvar(:,7))];
% bar(SLvarbaryoung)
% xticklabels({'Level','Incline','Decline'})
% legend('0.75 m/s','selfpace')
% title('SLvar')
% 
% figure
% SWvarbaryoung = [nanmean(aveSWvar(:,2)),nanmean(aveSWvar(:,5));nanmean(aveSWvar(:,8)),nanmean(aveSWvar(:,6));nanmean(aveSWvar(:,9)),nanmean(aveSWvar(:,7))];
% bar(SWvarbaryoung)
% xticklabels({'Level','Incline','Decline'})
% legend('0.75 m/s','selfpace')
% title('SWvar')
% 
% figure
% SWbaryoung = [nanmean(aveSW(:,2)),nanmean(aveSW(:,5));nanmean(aveSW(:,8)),nanmean(aveSW(:,6));nanmean(aveSW(:,9)),nanmean(aveSW(:,7))];
% bar(SWbaryoung)
% xticklabels({'Level','Incline','Decline'})
% legend('0.75 m/s','selfpace')
% title('SW')
% 
% figure
% SPbaryoung = [nanmean(GRP.st(:,2)),nanmean(GRP.st(:,5));nanmean(GRP.st(:,8)),nanmean(GRP.st(:,6));nanmean(GRP.st(:,9)),nanmean(GRP.st(:,7))];
% bar(SPbaryoung)
% xticklabels({'Level','Incline','Decline'})
% legend('0.75 m/s','selfpace')
% title('SP')
% 
% figure
% SPvarbaryoung = [nanmean(GRP.stvar(:,2)),nanmean(GRP.stvar(:,5));nanmean(GRP.stvar(:,8)),nanmean(GRP.stvar(:,6));nanmean(GRP.stvar(:,9)),nanmean(GRP.stvar(:,7))];
% bar(SPvarbaryoung)
% xticklabels({'Level','Incline','Decline'})
% legend('0.75 m/s','selfpace')
% title('SPvar')
