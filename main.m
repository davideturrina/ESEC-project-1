close all; 
clear;   
clc;

%% DATA EXTRACTION

% Demand data - demand in kWh 
demand = 'demand_data_hxh_8784h.csv';
demand = readtable(demand, 'PreserveVariableNames', true);
demand = demand(1:8784, 1:2);

% Wind turbine power curve data (Power in kW: second column)
wind_powercurve = xlsread('turbine_power_curve_5_MW.xlsx','Sheet1','B2:E32');
wind_powercurve= wind_powercurve(:,1:2);

% Solar and wind data for each hour in Orlando
solar_wind_hxh= 'solar_and_wind_data_hxh.csv'; 
solar_wind_hxh= readtable(solar_wind_hxh, 'PreserveVariableNames', true);
solar_wind_hxh= solar_wind_hxh(1:8784, 1:7);

% Solar efficiency and panel area
solar_efficiency = 0.197;
Area=2.80;

% Extract numerical values for solar and wind data
sol5_numeric = double(solar_wind_hxh{:,5});
wind7_numeric = double(solar_wind_hxh{:, 7});

%% VARIOUS PLOTS

% Plot demand data
figure;
plot(demand{:,1}, demand{:,2}, 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Demand (kW)');
title('Energy Demand');
grid on;

% Plot wind power curve
figure;
plot(wind_powercurve(:,1), wind_powercurve(:,2), 'LineWidth', 1.5);
xlabel('Wind Speed (m/s)');
ylabel('Power Output (kW)');
title('Wind Turbine Power Curve');
grid on;

% Plot solar and wind data with dual y-axes
time = 1:height(solar_wind_hxh); % Create time axis (hours)
figure;

% Y-axis for solar data
yyaxis left;
plot(time, solar_wind_hxh{:,5}, 'DisplayName', 'Solar Data', 'LineWidth', 1.5);
ylabel('Solar Power (kW)');
ylim([0 max(solar_wind_hxh{:,5})*1.1]); % Add 10% upper margin

% Y-axis for wind data
yyaxis right;
plot(time, solar_wind_hxh{:,7}, 'DisplayName', 'Wind Data', 'LineWidth', 1.5);
ylabel('Wind Power (kW)');
ylim([0 max(solar_wind_hxh{:,7})*1.1]); % Add 10% upper margin

% Titles and X-axis
xlabel('Time (hours)');
title('Solar and Wind Power Data');
legend('Solar Power', 'Wind Power');
grid on;

%% POWER ANALYSIS

% Calculate energy output for one solar panel
OnePanelEnergyOUTPUT = sol5_numeric .* Area .* solar_efficiency;

% Plot energy output for one solar panel
figure; % Create a new figure for the plot
plot(OnePanelEnergyOUTPUT);
xlabel('Index');
ylabel('Energy (E)');
title('Energy Output (E) Plot');
grid on;

% Calculate energy output for one wind turbine
% Input data: power table (power in kW) as a function of wind speed (in m/s)
wind_powercurve(:,1); % Wind speeds
wind_powercurve(:,2); % Corresponding power outputs

% Calculate energy produced for each hour
energy_per_hour = zeros(length(wind7_numeric), 1); % Preallocate energy vector

for i = 1:length(wind7_numeric)
    % Find the index of the closest wind speed in the power table
    %AFONSO: do you think that there is another way ? there is an
    %aproximation but I think this is the best and easiest way. Otherwise
    %we can do an interpolation of the power curve and..
    [~, idx] = min(abs(wind_powercurve(:,1) - wind7_numeric(i)));
    % Get the corresponding power output
    power = wind_powercurve(idx, 2); % Correctly access the power value
    % Energy produced in one hour (kWh)
    energy_per_hour(i) = power * 1; % Power (kW) * Time (1 hour)
end

% Total energy produced in a year
one_turbine_total_energy = sum(energy_per_hour);

% Display results
disp(['Energia totale annuale per una turbina: ', num2str(one_turbine_total_energy), ' kWh']);


%% W\ INTERPOL
% wind_speeds = wind_powercurve(:,1); 
% power_output = wind_powercurve(:,2); 
% interpolated_power = interp1(wind_speeds, power_output, wind7_numeric, 'linear', 'extrap');
% energy_per_hour = interpolated_power; 
% one_turbine_total_energy = sum(energy_per_hour);
% disp(['Energia totale prodotta in un anno con interpolazione: ', num2str(one_turbine_total_energy), ' kWh']);

