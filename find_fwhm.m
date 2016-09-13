function fwhm=find_fwhm(tuning_curve)

deg_per_ind=360/length(tuning_curve);
curve_max=max(tuning_curve);
curve_min=min(tuning_curve);

max_index=find(tuning_curve==curve_max);

half_max=curve_max-(curve_max-curve_min)/2;

epsilon=1;

while 1

    fwhm=epsilon*2*deg_per_ind;

    up_index=max_index+epsilon;
    down_index=max_index-epsilon;

    if up_index>length(tuning_curve);
        up_index=up_index-length(tuning_curve);
    end

    if down_index<1
        down_index=length(tuning_curve)-abs(down_index);
    end

    if mean([tuning_curve(up_index) tuning_curve(down_index)])<half_max
        break;
    end

    epsilon=epsilon+1;

end



