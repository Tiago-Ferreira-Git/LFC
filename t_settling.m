function time_settle = t_settling(t,t_transient,freq_lim,freq,h,t_shift)

    
    index = find(abs(freq) < freq_lim);
    
    mask = find(t >= t_transient);
    if length(mask) == 0
        mask;
    end
    mask = mask(1);
    %Ignore initial transient
    mask = index > mask;
    index(~mask) = [];
    %index = index(1):h:index(end);
    
    for j = 1:length(index)- round(5/h)
        if(index(j+round(5/h)) - index(j) == round(5/h))
            index = index(j);
            break
        end
    end
    
    if length(index) ~=1
        error 'Could not find time settling'
    end
    
    time_settle =  t(index)-t_shift;
    if time_settle > 3600
        time_settle = 0;
    end
        
end

