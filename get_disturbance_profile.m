function w = get_disturbance_profile(w,h,n_areas,simulation_hours)
    w_load_hour = load('data\w_load.mat');
    w_load_hour = w_load_hour.w_load(:,1:24);
    
    w_ren_hour = load('data\w_ren.mat');
    w_ren_hour = w_ren_hour.w_ren(:,1:24);

    if size(w_load_hour,1) < n_areas 
        w_load_hour = [w_load_hour; w_load_hour(1:n_areas-size(w_load_hour,1),:)];
        w_ren_hour = [w_ren_hour; w_ren_hour(1:n_areas-size(w_ren_hour,1),:)];
    else 
        w_load_hour = w_load_hour(1:n_areas,:);
        w_ren_hour = w_ren_hour(1:n_areas,:);
    end


    if 24 < simulation_hours 
        w_load_hour = repmat(w_load_hour,1,floor(simulation_hours/24)+1);
        w_ren_hour =  repmat(w_ren_hour,1,floor(simulation_hours/24)+1);
    end


    k_ = 1:1440:size(w,2);
    for hour=1:simulation_hours
        w(1:2:end,k_(hour):k_(hour+1)) = repmat(w_load_hour(:,hour),1,1441)*100;
        w(2:2:end,k_(hour):k_(hour+1)) = repmat(w_ren_hour(:,hour),1,1441);
    end
end

