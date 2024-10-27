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
    

    % Retrieve Load in MW
    y_meas = table{:,3};
    y_forecast = table{:,2};



    %last data point time
    last_date = table{end,1}{1};
    last_date = strsplit(last_date,'-');
    last_date = datetime(last_date(2),"InputFormat", "dd.MM.yyyy HH:mm");

    
    % Resample if needed / Obtain the difference array
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
        y_meas = y_meas(1:24*365+1)';
        
        y_forecast = y_forecast(1:24*365+1);
    end
    

    %Substituting the Nan Values by their preceding one 

    y_forecast = fillmissing(y_forecast,'previous');
    y_meas =  fillmissing(y_meas,'previous');
    
    normalization = y_forecast(1);

    y_forecast = y_forecast./normalization;
    y_meas = y_meas./normalization;

    w_forecast(j,:) = y_forecast;
    w_measured(j,:) = y_meas;

end
%%

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
ylabel('Load - Forecast')
hold off;
savefig('forecast.fig');
set(gcf,'renderer','Painters');
saveas(gca,'load_forecast.eps','epsc');


figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
stairs(t(:,mask),w_measured(:,mask)','LineWidth',1.5)
xlabel('Time (s)')
ylabel('Load  - Measured')
hold off;
savefig('measured.fig');
set(gcf,'renderer','Painters');
saveas(gca,'load_measured.eps','epsc');




figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
stairs(t(:,mask),w_measured(:,mask)' - w_forecast(:,mask)','LineWidth',1.5)
xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
xticks(0:3600:(24-1)*3600)
xticklabels(sprintfc('%d', 0:(24-1)))
xlim([0 (24-1)*3600])
ylabel('Normalized Load Difference')
savefig('difference.fig');
set(gcf,'renderer','Painters');
saveas(gca,'load_difference.eps','epsc');
saveas(gca,'load_difference.png','png');



load_.description = "This data is normalized! Hourly data for one year!";
load_.measured = w_measured;
load_.forecast = w_forecast;



save('../../load_profile.mat','load_')
