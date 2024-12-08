close all; 
clear;   
clc;

%% DATA EXTRACTION

population = 2273800;

demandxcapita = readtable('demand_data_hxh_8784h.csv', 'PreserveVariableNames', true);
demandxcapita = demandxcapita(1:8784, 1:2);
demand_total = array2table(demandxcapita{:,1:2}, 'VariableNames', {'Hour', 'TotalDemand'});

demand_total.TotalDemand = demand_total.TotalDemand * population;



% Wind turbine power curve data (Power in kW: second column)
wind_powercurve = xlsread('turbine_power_curve_5_MW.xlsx', 'Sheet1', 'B2:E32');
wind_powercurve = wind_powercurve(:, 1:2);

% Solar and wind data for each hour in Orlando
solar_wind_data = readtable('solar_and_wind_data_hxh.csv', 'PreserveVariableNames', true);
solar_wind_data = double(table2array(solar_wind_data(1:8784, 1:7)));
    
% Solar efficiency and panel area
solar_efficiency = 0.197*0.96;              %cell+inverter
panel_area = 2.80; %m2
solar_data5 = solar_wind_data(:, 5) / 1000;                 % Solar irradiance was in Wh, now is in kWh

solar_data5(solar_wind_data(:,6) > 80) = 0;                 %accounnting for shadows beyond 90-10 degrees


% Extract numerical values for solar and wind data
wind_data7 = solar_wind_data(:,7);
%accounting for roughness
h1 = 50; h2 = 140; z0 = 1.6;
wind_data7 = wind_data7 .* (log(h2 / z0) / log(h1 / z0)); 


%% Energy output single panel
solar_output_matrix = solar_data5 .* panel_area .* solar_efficiency;
total_solar_output = sum(solar_output_matrix);
disp(['Total energy for one PV (2.80m²) in a year: ', num2str(total_solar_output), ' kWh']);

%% Energy output single turbine (with interpolation )
 wind_speeds = wind_powercurve(:,1); 
 power_output = wind_powercurve(:,2); 
 interpolated_power = interp1(wind_speeds, power_output, wind_data7, 'linear', 'extrap');
 energy_per_hour = interpolated_power; 
 one_turbine_total_energy = sum(energy_per_hour);
 disp(['Total energy produce for one turbine (w\ interpol): ', num2str(one_turbine_total_energy), ' kWh']);


%% PLOTS

% Plot della domanda totale
figure;
plot(demand_total.Hour, demand_total.TotalDemand, 'LineWidth', 1.5);
xlabel('Time (Hours)');
ylabel('Total Demand (kW)');
title('Energy Demand Over Time');
grid on;

% Plot wind power curve
figure;
plot(wind_powercurve(:,1), wind_powercurve(:,2), 'LineWidth', 1.5);
xlabel('Wind Speed (m/s)');
ylabel('Power Output (kW)');
title('Wind Turbine Power Curve');
grid on;

% Plot solar and wind data with dual y-axes
time = 1:height(solar_wind_data); 
figure;
yyaxis left;
plot(time, solar_wind_data(:,5), 'DisplayName', 'Solar Data', 'LineWidth', 1.5);
ylabel('Solar Power (kW)');
ylim([0 max(solar_wind_data(:,5))*1.1]); 
yyaxis right;
plot(time, solar_wind_data(:,7), 'DisplayName', 'Wind Data', 'LineWidth', 1.5);
ylabel('Wind Power (kW)');
ylim([0 max(solar_wind_data(:,7))*1.1]); 
xlabel('Time (hours)');
title('Solar and Wind Power Data');
legend('Solar Power', 'Wind Power');
grid on;







