function [STV_OP1,STV_OP1_median,STV_OP1_steadypoint_all,STV_OP1_steadypoint] = stdmethod(stridespeed)


    
STV_OP1 = movstd(stridespeed(:),[0 5]);

STV_OP1_median = median(STV_OP1)-(0.25*median(STV_OP1));
 
% STV_OP1_median = 0.006; % for fixed speed 


STV_OP1_steadypoint_first = find(STV_OP1<STV_OP1_median);

STV_OP1_steadypoint = find(STV_OP1<STV_OP1_median,1);


for i= 1:5
   STV_OP1_steadypoint_first(:,i+1) = STV_OP1_steadypoint_first(:,1)+ i;   
end

STV_OP1_steadypoint_all = unique(STV_OP1_steadypoint_first);
% STV_OP1_steadypoint_all = STV_OP1_steadypoint_all';

indices = find(STV_OP1_steadypoint_all>=394);
STV_OP1_steadypoint_all(indices) = []; 
end 

    
