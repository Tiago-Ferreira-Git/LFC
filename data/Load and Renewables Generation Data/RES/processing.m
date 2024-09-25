clear all;
close all;


myDir = pwd; %gets directory
myFiles = dir(fullfile(myDir,'2022/','*.csv')); %gets all wav files in struct

w_ren = zeros(length(myFiles),24*365+1);
w_ren_profile = w_ren;
teste = w_ren;


for j = 1:length(myFiles)

    baseFileName = fullfile(myDir,'2022/',myFiles(j).name);
    
    table = readtable(baseFileName);


    % Retrive sampling time

    string = table{1,1}{1};
    string = strsplit(string,'-');
    
    first_date = datetime(string(1),"InputFormat", "dd.MM.yyyy HH:mm");


    time_elapsed = datetime(string(2),"InputFormat", "dd.MM.yyyy HH:mm") - datetime(string(1),"InputFormat", "dd.MM.yyyy HH:mm");
    ts = seconds(time_elapsed);

    fs = 1/ts;
    

    % Retrieve Solar in MW
    y_meas = table{:,3};
    y_forecast = table{:,2};

    if ~isa(y_meas,'double') || ~isa(y_forecast,'double') || all(isnan(y_forecast)) || all(isnan(y_meas))
        % Retrieve Wind and use it instead in MW
        y_meas = table{:,7};
        y_forecast = table{:,6};
        if ~isa(y_meas,'double') || ~isa(y_forecast,'double') || all(isnan(y_forecast)) || all(isnan(y_meas))
            continue
        end
    end



    %last data point time
    last_date = table{end,1}{1};
    last_date = strsplit(last_date,'-');
    last_date = datetime(last_date(2),"InputFormat", "dd.MM.yyyy HH:mm");

    
    % Resample if needed / Obtain the difference array
    if ts ~= 3600
        y_meas = resample(y_meas,ts,3600);
        y_meas = y_meas(1:24*365+1);


        y_forecast = resample(y_forecast,ts,3600);
        y_forecast = y_forecast(1:24*365+1);

        w_ren(j,:) = y_forecast;

        teste(j,:) = y_meas;
        
        % y = filter([-1 1],[1 0],y_resampled);
        % y(1) = 0;

    else
        y_meas = y_meas(1:24*365+1)';
        y_meas = y_meas(1:24*365+1)';
        
        y_forecast = y_forecast(1:24*365+1);

        w_ren(j,:) = y_forecast;
        teste(j,:) = y_meas;
    end
    

    %Substituting the Nan Values by their preceding one 


    while(any(isnan(y_meas)))
         mask = circshift(isnan(y_meas),-1);
         y_meas(isnan(y_meas)) = y_meas(mask);
    end


    while(any(isnan(y_forecast)))            
        mask = circshift(isnan(y_forecast),-1);
        y_forecast(isnan(y_forecast)) = y_forecast(mask);
    end


    w_ren_profile(j,:) = 1 + (y_meas-y_forecast)'/mean(y_forecast);


    time_elapsed = size(w_ren,2);

    t = 0:3600:(time_elapsed-1)*3600;
    mask = t < 3600*24;
        
    figure
    plot(t,w_ren(j,:))
    xlabel('Time (s)')
    ylabel('Solar - pu - Forecast')
    
    
    figure
    plot(t,teste(j,:))
    xlabel('Time (s)')
    ylabel('Soalr - MW - Measured')

end



%%

figure
stairs(t(mask),w_ren_profile(:,mask)')
xlabel('Time (s)')
ylabel('Load - $1 + \frac{measured - forecasted}{mean(forecasted)}$','Interpreter','latex')


save('../../ren_profile.mat','w_ren_profile')