function [time_n, data_n, data_nft] = normalize_gaitcycle2(time, data, GEidx, npts)

% data and GEidx must be a matrix, not a structure

time_n = nan(npts,length(GEidx));

if size(data,2) == 1
    data_n = nan(npts, length(GEidx));
    data_nft = nan(npts, length(GEidx));
else
    data_n = nan(npts, size(data,2), length(GEidx));
    data_nft = nan(npts, size(data,2), length(GEidx));
end

for sn = 1:length(GEidx)-1 %loop through strides
    
    % create linearly spaced time vector with npts
    time_n(:,sn) = linspace(time(GEidx(sn,1)),time(GEidx(sn,5)),npts)';
    %     temptime = time(GEidx(sn,1):GEidx(sn+1,1));
    %     temptime1 = temptime([1:2:100 101:length(temptime)]);
    
    if size(data,2) == 1
        data_n(:,sn) = interp1(time(GEidx(sn,1):GEidx(sn,5)), data(GEidx(sn,1):GEidx(sn,5),:), time_n(:,sn),'spline');
        %     data_temp(:,sn) = interp1(time(GEidx(sn,1):GEidx(sn+1,1)), data(GEidx(sn,1):GEidx(sn+1,1),:), temptime1,'spline');
        data_nft(:,sn) = interpft(data(GEidx(sn,1):GEidx(sn,5),:), npts);
    else
        % interp
        data_n(:,:,sn) = interp1(time(GEidx(sn,1):GEidx(sn,5)), data(GEidx(sn,1):GEidx(sn,5),:), time_n(:,sn),'spline');
        %     data_temp(:,:,sn) = interp1(time(GEidx(sn,1):GEidx(sn+1,1)), data(GEidx(sn,1):GEidx(sn+1,1),:), temptime1,'spline');
        data_nft(:,:,sn) = interpft(data(GEidx(sn,1):GEidx(sn,5),:), npts);
    end
end
