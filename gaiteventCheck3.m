function goodstrides = gaiteventCheck3(RHS, LTO, LHS, RTO, vRGRF, vLGRF)

colorleg = {'b' 'r'}; % blue for left, red for right

flag_plot = 0;
%%
if flag_plot
    figure
    plot(vRGRF,'r');
    hold on
    plot(vLGRF,'b');
    plot(RHS, zeros(size(RHS)), 'rx', LTO, zeros(size(LTO)), 'bo', LHS, zeros(size(LHS)), 'bx', RTO, zeros(size(RTO)), 'ro')
end
medianST = mean([median(diff(RHS)) median(diff(LTO)) median(diff(LHS)) median(diff(RTO))]);

%% check leg
evsorted = sort([RHS; LTO; LHS; RTO]);

for e = 1:length(evsorted)
    if vRGRF(evsorted(e)) < vLGRF(evsorted(e))
        leg(e) = 2; % right = 2
    else
        leg(e) = 1; % left = 1
    end
end

if flag_plot
    R = evsorted(leg == 2);
    L = evsorted(leg == 1);
    
    plot(R, zeros(size(R)), [colorleg{2} '.'], 'markersize',16)
    plot(L, zeros(size(L)), [colorleg{1} '.'], 'markersize',16)
end

%%

for i = 1:length(RHS)
    rhsi(i) = find(evsorted == RHS(i),1,'first');
end

for i = 1:length(LTO)
    ltoi(i) = find(evsorted == LTO(i),1,'first');
end

for i = 1:length(LHS)
    lhsi(i) = find(evsorted == LHS(i),1,'first');
end

for i = 1:length(RTO)
    rtoi(i) = find(evsorted == RTO(i),1,'first');
end

evtype = nan([length(evsorted) 1]);
evtype(rhsi,1) = 1;
evtype(ltoi,1) = 2;
evtype(lhsi,1) = 1;
evtype(rtoi,1) = 2;

%%
tempRHS = []; tempLTO = []; tempLHS = []; tempRTO = [];
for e = 1:length(evsorted)
    if evtype(e) == 1 && leg(e) == 2
        tempRHS = [tempRHS; evsorted(e)];
    elseif evtype(e) == 2 && leg(e) == 1
        tempLTO = [tempLTO; evsorted(e)];
    elseif evtype(e) == 1 && leg(e) == 1
        tempLHS = [tempLHS; evsorted(e)];
    elseif evtype(e) == 2 && leg(e) == 2
        tempRTO = [tempRTO; evsorted(e)];
    end 
end
%%
if flag_plot
    figure
    plot(vRGRF,'r');
    hold on
    
    plot(vLGRF,'b');
    plot(tempRHS, zeros(size(tempRHS)), 'rx', tempLTO, zeros(size(tempLTO)), 'bo', tempLHS, zeros(size(tempLHS)), 'bx', tempRTO, zeros(size(tempRTO)), 'ro')
end
%%
% 
% rhsi = 1; ltoi = 1; lhsi=1; rtoi=1;
% rhsct = 1; ltoct = 1; lhsct=1; rtoct=1;
% for e = 1:length(evsorted)
%     evsorted(e)
%     RHS(rhsi)
%     if evsorted(e) == RHS(rhsi) && leg(e) == 2
%         goodRHS(rhsct) = evsorted(e);
%         rhsi = rhsi+1;
%         rhsct = rhsct+1;
%     elseif evsorted(e) == RHS(rhsi) && leg(e) == 1
%         goodLHS(lhsct) = evsorted(e);
%         lhsi = lhsi+1;
%         lhsct = lhsct+1;
%     end
%     
%     evsorted(e)
%     LTO(ltoi)
%     if evsorted(e) == LTO(ltoi) && leg(e) == 1
%         goodLTO(ltoct) = evsorted(e);
%         ltoi = ltoi+1;
%         ltoct = ltoct+1;
%     elseif evsorted(e) == LTO(ltoi) && leg(e) == 2
%         goodRTO(rtoct) = evsorted(e);
%         rtoi = rtoi+1;
%         rtoct = rtoct+1;
%     end
%     
%     evsorted(e)
%     LHS(lhsi)
%     if evsorted(e) == LHS(lhsi) && leg(e) == 1
%         goodLHS(lhsct) = evsorted(e);
%         lhsi = lhsi+1;
%         lhsct = lhsct+1;
%     elseif evsorted(e) == LHS(lhsi) && leg(e) == 2
%         goodRHS(rhsct) = evsorted(e);
%         rhsi = rhsi+1;
%         rhsct = rhsct+1;
%     end
%     
%     evsorted(e)
%     RTO(rtoi)
%     if evsorted(e) == RTO(rtoi) && leg(e) == 2
%         goodRTO(rtoct) = evsorted(e);
%         rtoi = rtoi+1;
%         rtoct = rtoct+1;
%     elseif evsorted(e) == RTO(rtoi) && leg(e) == 1
%         goodLTO(ltoct) = evsorted(e);
%         ltoi = ltoi+1;
%         ltoct = ltoct+1;
%     end
% %     pause
%     disp('here')
% end
% 
% keyboard
%%
% if flag_plot
%     figure
%     plot(vRGRF,'r');
%     hold on
%     
%     plot(vLGRF,'b');
%     plot(goodRHS, zeros(size(goodRHS)), 'rx', goodLTO, zeros(size(goodRHS)), 'bo', goodLHS, zeros(size(LHSgoodLHStemp)), 'bx', goodRTO, zeros(size(goodRTO)), 'ro')
% end

%%
% rhsct = 1; lhsct = 1;
% rtoct = 1; ltoct = 1;
% for e = 1:length(evsorted)
%     if evtype(e) == 1 && leg(e) == 2
%         RHStemp(rhsct) = evsorted(e);
%         rhsct = rhsct + 1;
%     elseif evtype(e) == 1 && leg(e) == 1
%         LHStemp(lhsct) = evsorted(e);
%         lhsct = lhsct + 1;
%     elseif evtype(e) == 2 && leg(e) == 2
%         RTOtemp(rtoct) = evsorted(e); 
%         rtoct = rtoct + 1;
%     elseif evtype(e) == 2 && leg(e) == 1
%         LTOtemp(ltoct) = evsorted(e); 
%         ltoct = ltoct + 1;
%     end
% end

% RHSi = intersect(find(evtype == 1), find(leg == 2));
% LTOi = intersect(find(evtype == 2), find(leg == 1));
% LHSi = intersect(find(evtype == 1), find(leg == 1));
% RTOi = intersect(find(evtype == 2), find(leg == 2));
% 
% RHStemp = evsorted(RHSi);
% LTOtemp = evsorted(LTOi);
% LHStemp = evsorted(LHSi);
% RTOtemp = evsorted(RTOi);
% 
% if flag_plot
%     figure
%     plot(vRGRF,'r');
%     hold on
%     
%     plot(vLGRF,'b');
%     plot(RHStemp, zeros(size(RHStemp)), 'rx', LTOtemp, zeros(size(LTOtemp)), 'bo', LHStemp, zeros(size(LHStemp)), 'bx', RTOtemp, zeros(size(RTOtemp)), 'ro')
% end

% keyboard
%% Need to get into RHS LTO LHS RTO


% GEchar = [leg' evtype'];

ct=1; e=1;

% numstride = min([length(tempRHS) length(tempLTO) length(tempLHS) length(tempRTO)]);

rhsct = 1; lhsct = 1;
rtoct = 1; ltoct = 1;
rhs1ct = 1;
goodRHS1 = 0;
go=1;

while go
    % first column is leg, R = 2, L = 1
    % second column is HS = 1, TO = 2
    lookingforRHS = 1; lookingforLTO = 1; lookingforLHS = 1; lookingforRTO = 1;
    lookingforRHS1 = 1;
    
    while lookingforRHS && e <= length(evsorted)
        if e==1 && evtype(e) == 1 && leg(e) == 2 %RHS
            goodRHS = evsorted(e);
            %                 rhsct = rhsct+1;
            lookingforRHS = 0;
        else 
            e=e+1;
        end
        
        if e>1 && e <= length(evsorted) %RHS
            if evtype(e) == 1 && leg(e) == 2 && evsorted(e) < goodRHS1
                e=e+1;
            else
                goodRHS = goodRHS1;
%                 rhsct = rhsct+1;
                lookingforRHS = 0;
            end
        end
    end

    while lookingforLTO && e <= length(evsorted)
        if evtype(e) == 2 && leg(e) == 1 %LTO
            goodLTO = evsorted(e);
%             ltoct = ltoct+1;
            lookingforLTO = 0;
        else
            e=e+1;
        end
    end
    
    
    while lookingforLHS && e <= length(evsorted)
        if evtype(e) == 1 && leg(e) == 1 %LHS
            goodLHS = evsorted(e);
%             lhsct = lhsct+1;
            lookingforLHS = 0;
        else
            e=e+1;
        end
    end
    
    while lookingforRTO && e <= length(evsorted)
        if evtype(e) == 2 && leg(e) == 2 %RTO
            goodRTO = evsorted(e);
%             rtoct = rtoct+1;
            lookingforRTO = 0;
        else
            e=e+1;
        end
    end
    
    while lookingforRHS1 && e <= length(evsorted)
        if evtype(e) == 1 && leg(e) == 2 %RHS
            goodRHS1 = evsorted(e);
%             rhs1ct = rhs1ct+1;
            lookingforRHS1 = 0;
        else
            e=e+1;
        end
    end
    
    tempstrides(ct,:) = [goodRHS goodLTO goodLHS goodRTO goodRHS1];
    ct = ct+1;
    
    if e >= length(evsorted), go = 0; end
    
end


%%
% tempstrides = [RHStemp LTOtemp LHStemp RTOtemp];
if flag_plot
    figure
    plot(vRGRF,'r');
    hold on
    
    plot(vLGRF,'b');
    plot(tempstrides(:,1), zeros(size(tempstrides(:,1))), 'rx', tempstrides(:,2), zeros(size(tempstrides(:,1))), 'bo', tempstrides(:,3), zeros(size(tempstrides(:,1))), 'bx', tempstrides(:,4), zeros(size(tempstrides(:,1))), 'ro')
end
%%
% check vLGRF reached < 25 N between LTO and LHS
% if not, do not include
ct=1;
for sn = 1:(length(tempstrides)-1)
    tempRHS0 = tempstrides(sn,1);
    tempLTO = tempstrides(sn,2);
    tempLHS = tempstrides(sn,3);
    tempRTO = tempstrides(sn,4);
    tempRHS1 = tempstrides(sn,5);
    if mean(abs(vLGRF(tempLTO:tempLHS))) < 25 && mean(abs(vRGRF(tempRTO:tempRHS1))) < 25
        tempstrides1(ct,:) = tempstrides(sn,:);
        ct=ct+1;
    end
end

%%
if flag_plot
    figure
    plot(vRGRF,'r');
    hold on
    
    plot(vLGRF,'b');
    plot(tempstrides1(:,1), zeros(size(tempstrides1(:,1))), 'rx', tempstrides1(:,2), zeros(size(tempstrides1(:,1))), 'bo', tempstrides1(:,3), zeros(size(tempstrides1(:,1))), 'bx', tempstrides1(:,4), zeros(size(tempstrides1(:,1))), 'ro',tempstrides1(:,5), zeros(size(tempstrides1(:,5))), 'rx')
end
%%
% check stride times are reasonable, ~0.6*medianST to 1.4*medianST
for sn = 1:(length(tempstrides1)-1)
    if tempstrides1(sn+1,1)-tempstrides1(sn,1) < 0.7*medianST || tempstrides1(sn+1,1)-tempstrides1(sn,1) > 1.4*medianST
        tempstrides1(sn,2:4) = NaN;
    end
end
%%
if flag_plot
    figure
    plot(vRGRF,'r');
    hold on
    
    plot(vLGRF,'b');
    plot(tempstrides1(:,1), zeros(size(tempstrides1(:,1))), 'rx', tempstrides1(:,2), zeros(size(tempstrides1(:,1))), 'bo', tempstrides1(:,3), zeros(size(tempstrides1(:,1))), 'bx', tempstrides1(:,4), zeros(size(tempstrides1(:,1))), 'ro')
end
%%
% check no singular RHS
ct=1;
for sn = 1:(length(tempstrides1)-1)
    if sum(isnan([tempstrides1(sn,:) tempstrides1(sn+1,1)])) == 0
        goodstrides(ct,:) = [tempstrides1(sn,:)];
        ct=ct+1;
    end
end

%%
if flag_plot
    figure
    plot(vRGRF,'r');
    hold on
    
    plot(vLGRF,'b');
    plot(goodstrides(:,1), zeros(size(goodstrides(:,1))), 'rx', goodstrides(:,2), zeros(size(goodstrides(:,1))), 'bo', goodstrides(:,3), zeros(size(goodstrides(:,1))), 'bx', goodstrides(:,4), zeros(size(goodstrides(:,1))), 'ro', goodstrides(:,5), zeros(size(goodstrides(:,1))), 'rx')
end

