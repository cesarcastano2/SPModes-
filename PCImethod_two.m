function [PCI_percent,step_stride_norm,PCI_CV_OP1,PCI_CV_OP2,ROI] = PCImethod_two(steptime,stridetime)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
step_stride_norm = 360.*(steptime(:,1)./stridetime(:,2));

for i = 1:length(step_stride_norm)-1
PCI_percent(i,1) = (100*(std(step_stride_norm(1:i+1))/mean(step_stride_norm(1:i+1))))+(100*(mean(abs(step_stride_norm(1:i+1)-180)/180)));
end 
PCI_percent_flipped = flipud(PCI_percent);
% Option 1
for i = 1:length(PCI_percent_flipped)-15
PCI_CV_OP1(i) = std(PCI_percent_flipped(i:(i+15)))./mean(PCI_percent_flipped(i:(i+15)));
end
PCI_CV_OP1_size = length(stridetime)-(find(PCI_CV_OP1>0.05));
% Option 2
for i = 14:length(PCI_percent_flipped)-1
PCI_CV_OP2(:) =movstd(PCI_percent_flipped,[0 (i)])./movmean(PCI_percent_flipped,[0 (i)]);
if find(PCI_CV_OP2(1:(i+1))>0.05)==1
    PCI_CV_OP2_size(:) = length(stridetime)-length(PCI_CV_OP2);
else 
    PCI_CV_OP2_size(1) = 0;
end
% end
PCI_CV_OP2_size(1) = 0;
% Highest value for region of interest 
if PCI_CV_OP1_size(1)>= PCI_CV_OP2_size(1)
   ROI(1) = PCI_CV_OP1_size(1);
else 
    ROI(1) = PCI_CV_OP2_size(1);
end

% POS = length(stridetime)-ROI;
% 
% for i = 14:-2:2
%     POS_finetuning = movstd(PCI_percent_flipped((i-POS):POS,[0 (i)]))./movmean(PCI_percent_flipped((i-POS):POS,[0 (i)]));
% end



% PCI_CV_steadystate = (length(find(PCI_CV>0.05)))+1;
end