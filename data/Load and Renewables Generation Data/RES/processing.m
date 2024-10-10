clear all;
close all;


myDir = pwd; %gets directory
myFiles = dir(fullfile(myDir,'2022/','*.csv')); %gets all wav files in struct

w_forecast = zeros(length(myFiles),24*365+1);
w_measured = w_forecast;


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

    
    if ts ~= 3600

        if ts > 3600
            y_meas = resample(y_meas,ts,3600);
            y_meas = y_meas(1:24*365+1);
    
    
            y_forecast = resample(y_forecast,ts,3600);
            y_forecast = y_forecast(1:24*365+1);
        else
            y_forecast = y_forecast(1:3600/ts:end);
            y_meas = y_meas(1:3600/ts:end);
        
        end

    else
        y_meas = y_meas(1:24*365+1)';
        y_forecast = y_forecast(1:24*365+1)';
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
    
    
    %normalization in Renewables will be its maximum value in the year,
    %such that the user can specify it without needing to change the whole
    %year profile

    
    normalization = max(y_forecast);



    y_forecast = y_forecast./normalization;
    y_meas = y_meas./normalization;

    w_forecast(j,:) = y_forecast;
    w_measured(j,:) = y_meas;


end

time_elapsed = size(w_forecast,2);

t = 0:3600:(time_elapsed-1)*3600;
mask = t < 3600*24*1;

figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
stairs(t(:,mask),w_forecast(:,mask)','LineWidth',1.5)
xlabel('Time (s)')
ylabel('Res - Forecast')
hold off;
savefig('forecast.fig');
set(gcf,'renderer','Painters');
saveas(gca,'forecast.eps','epsc');


figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
stairs(t(:,mask),w_measured(:,mask)','LineWidth',1.5)
xlabel('Time (s)')
ylabel('Res  - Measured')
hold off;
savefig('measured.fig');
set(gcf,'renderer','Painters');
saveas(gca,'measured.eps','epsc');




figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
stairs(t(:,mask),w_measured(:,mask)' - w_forecast(:,mask)','LineWidth',1.5)
xlabel('Time (s)')
ylabel('Res  - Difference')
savefig('difference.fig');
set(gcf,'renderer','Painters');
saveas(gca,'difference.eps','epsc');



res.description = "This data is normalized! Hourly data for one year!";
res.measured = w_measured;
res.forecast = w_forecast;





save('../../res_profile.mat','res')
