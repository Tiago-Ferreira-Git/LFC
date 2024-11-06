function time_settle = t_settling(t,t_transient,freq_lim,freq,h,t_shift)
%t_settling   A function to determine the time that frequency takes to be regulated to omega_lim .
%
%   Inputs:
%       t           - Time vector (in seconds).
%       t_transient - Time (in seconds) to ignore at the beginning of the simulation 
%                     to allow transient effects to settle.
%       freq_lim    - Frequency deviation limit, specifying the tolerance for 
%                     considering the system as settled.
%       freq        - Vector of frequency values over time, from which settling time 
%                     is calculated.
%       h           - Time step in seconds for the frequency data.
%       t_shift     - Time shift applied to ensure the result falls within the 
%                     desired time range [0, 3600] seconds.
%
%   Output:
%       time_settle - Settling time (in seconds) when the frequency deviation first 
%                     remains within freq_lim for a consecutive period of time, 
%                     adjusted by t_shift. Returns 0 if time settling is not achieved 
%                     within the simulation period.
    
    index = find(abs(freq) < freq_lim);
    
    mask = find(t >= t_transient);
    if isempty(mask)
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

