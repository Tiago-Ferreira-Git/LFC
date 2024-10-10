%The IEEE118 has only forecast data that's why we need to have forecast
%mesaurements. You have to run processing before this 


clear all;
close all;

A = readtable('Generators.csv');

idx = (strcmp(A.Category,'Solar'));

A = A(idx,:);

A.MaxCapacity_MW_ = str2double(A.MaxCapacity_MW_);

bus = A.NodeOfConnection;

bus = char(bus);
bus = bus(:,5:end);
bus = str2num(bus);

max(A.MaxCapacity_MW_)

buses = sort(unique(bus));

for i = 1:length(buses)
    mask = ismember(bus,buses(i));
    bus(mask) = i;
end


myDir = pwd; %gets directory
forecast_files = dir(fullfile(myDir,'Forecast/','*.csv')); 
measured_files = dir(fullfile(myDir,'Measured/','*.csv')); 

w_forecast = zeros(length(buses),24*365+1);
w_measured = w_forecast;




for j = 1:length(forecast_files)


    %Get forecast
    baseFileName = fullfile(myDir,'Forecast/',forecast_files(j).name);
    table = readtable(baseFileName);

    
    first_date = datetime(table{1,1},"InputFormat", "dd.MM.yyyy HH:mm");
    time_elapsed = datetime(table{2,1},"InputFormat", "dd.MM.yyyy HH:mm") - first_date;
    ts = seconds(time_elapsed);

    fs = 1/ts;
    

    % Retrieve Solar in MW
    y_forecast = table{:,2};



    %Get forecast
    baseFileName = fullfile(myDir,'Measured/',measured_files(j).name);
    table = readtable(baseFileName);

    
    first_date = datetime(table{1,1},"InputFormat", "dd.MM.yyyy HH:mm");
    time_elapsed = datetime(table{2,1},"InputFormat", "dd.MM.yyyy HH:mm") - first_date;
    ts = seconds(time_elapsed);

    fs = 1/ts;
    

    % Retrieve Solar in MW
    y_meas = table{:,2};

    % figure
    % hold on
    % plot(y_forecast);
    % plot(y_meas);
    % legend({'Forecast','Measure'},'Location','best');

    y_meas = y_meas(1:24*365+1);
    y_forecast = y_forecast(1:24*365+1);

    w_forecast(bus(j),:) = w_forecast(bus(j),:) + y_forecast';
    w_measured(bus(j),:) = w_measured(bus(j),:) + y_meas';
    j;
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
% savefig('forecast.fig');
% set(gcf,'renderer','Painters');
% saveas(gca,'forecast.eps','epsc');


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
% savefig('measured.fig');
% set(gcf,'renderer','Painters');
% saveas(gca,'measured.eps','epsc');

%%

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

res.bus = buses;


save('../../../res_profile_118.mat','res')

