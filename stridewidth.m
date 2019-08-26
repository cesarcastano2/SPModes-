function stridewidth_leg = stridewidth(x_lead, x_trail, GEidx1, GEidx2)


for i = 1:length(GEidx1)-1
stridewidth_leg(i) = (x_lead(GEidx1(i))+x_lead(GEidx1(i+1))/2) - x_trail(GEidx2(i)) ;
end 
