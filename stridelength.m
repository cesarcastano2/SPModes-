function stridelength_leg = stridelength(y_lead, y_trail, GEidx1, GEidx2,Speed, Time, pitch)

% for i = 1:length(Speed)
% sl(i) = ((y_lead(GEidx1(i)) - y_trail(GEidx2(i))) + (Speed(i)*Time(i)))/cosd(pitch) ;
% end

stridelength_leg = ((y_lead(GEidx1) - y_trail(GEidx2))/cosd(pitch)) + (Speed.*Time) ;
end 
