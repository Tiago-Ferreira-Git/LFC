clear all;
close all;


myDir = pwd; %gets directory
myFiles = dir(fullfile(myDir,'2022/','*.csv')); %gets all wav files in struct

w_load = zeros(length(myFiles),24*365+1);
w_load_profile = w_load;
teste = w_load;


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
    

    % Retrieve Load in MW
    y_meas = table{:,3};
    y_forecast = table{:,2};



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

        w_load(j,:) = y_forecast;

        teste(j,:) = y_meas;
        
        % y = filter([-1 1],[1 0],y_resampled);
        % y(1) = 0;

    else
        y_meas = y_meas(1:24*365+1)';
        y_meas = y_meas(1:24*365+1)';
        
        y_forecast = y_forecast(1:24*365+1);

        w_load(j,:) = y_forecast;
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


    w_load_profile(j,:) = y_meas'/mean(y_meas);
    %w_load_profile(j,:) = y_forecast'/y_forecast(1);

end
%%

time_elapsed = size(w_load,2);

t = 0:3600:(time_elapsed-1)*3600;
mask = t < 3600*24*1;

figure
plot(t,w_load)
xlabel('Time (s)')
ylabel('Load - pu - Forecast')


figure
plot(t,teste)
xlabel('Time (s)')
ylabel('Load - MW - Measured')

%%


figure('Position',4*[0 0 2*192 144]); % Nice aspect ratio for double column
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex') % Latex style axis

stairs(t(mask),w_load_profile(:,mask)','LineWidth',1.5)
xlabel('Time (s)')
ylabel('Load - $\frac{measured}{mean(measured)}$','Interpreter','latex')
hold off;
% Save figure to .fig and .eps formats
savefig('Load.fig');
set(gcf,'renderer','Painters');
saveas(gca,'Load.svg','svg');





save('../../load_profile.mat','w_load_profile')
