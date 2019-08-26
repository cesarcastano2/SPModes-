function [speed_strides] = compile_speed(speed,s,c)

     if c == 1
     speed_strides.level_050(s,:) = speed(:);
     end
     if c == 2
     speed_strides.level_100(s,:) = speed(:);
     end
     if c == 3
     speed_strides.level_150(s,:) = speed(:);
     end
     
     if c == 4
     speed_strides.incline_050(s,:) =  speed(:);
     end
     if c == 5
     speed_strides.incline_100(s,:) =  speed(:);
     end
     if c == 6
     speed_strides.incline_150(s,:) =  speed(:);
     end
     
     if c == 7
     speed_strides.decline_050(s,:) =  speed(:);
     end
     if c == 8
     speed_strides.decline_100(s,:) =  speed(:);
     end
     if c == 9
     speed_strides.decline_150(s,:) =  speed(:);
     end
end

