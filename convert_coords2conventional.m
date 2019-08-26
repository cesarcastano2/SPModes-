function [converted_data] = convert_coords2conventional(datain)

% conventional coords
% +x = right, -x = left
% +y = fowards, -y = backwards
% +z = up, -z = down

% axes for Motek
% +x = right, -x = left
% +y = up, -y = down
% +z = backwards, -z = forwards

converted_data = [datain(:,1) -datain(:,3) datain(:,2)];
