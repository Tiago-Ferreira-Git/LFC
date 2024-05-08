function w = get_disturbance_profile(w,h,n_areas)
    w_load_hour = load('data\w_load.mat');
    w_load_hour = w_load_hour.w_load;

    w_ren_hour = load('data\w_ren.mat');
    w_ren_hour = w_ren_hour.w_ren;

    if size(w_load_hour,1) < n_areas 
        w_load_hour = [w_load_hour; w_load_hour(1:n_areas-size(w_load_hour,1),:)];
        w_ren_hour = [w_ren_hour; w_ren_hour(1:n_areas-size(w_ren_hour,1),:)];
    else 
        w_load_hour = w_load_hour(1:n_areas,:);
        w_ren_hour = w_ren_hour(1:n_areas,:);
    end


    hour = 1;
    i = 0;
    for k=1:size(w,2)
        t = k*h;
        if(hour*3600 < t)
            hour = hour + 1;
        end
        if(24 <= hour-i*24 )
            i = i+1;
            hour-i*24+1;
        end
        w(1:2:end,k) = w_load_hour(:,hour-i*24+1)*100;
        w(2:2:end,k) = w_ren_hour(:,hour-i*24+1);
    end
end

