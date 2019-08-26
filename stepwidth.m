function sw = stepwidth(x_lead, x_trail, GEidx)

sw0 = x_lead(GEidx) - x_trail(GEidx);
sw = abs(sw0); % because left heel strikes will be negative