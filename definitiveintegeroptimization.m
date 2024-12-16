close all;
clear;
clc;

%% INITIALIZATION
% Initialize variables to store optimal solution
storage_capacity = 0;
n_pv = 0;
n_wind = 0;
n_nuclear = 0;

%% DATA EXTRACTION
% Load population and area data
population = 2273800;
area = 11246.8; % km²
max_investment = population * 3000;

% Load demand data
% Each row represents hourly demand data for a single individual
% The second column is the demand per capita

demand_per_capita = readtable('demand_data_hxh_8784h.csv', 'PreserveVariableNames', true);
demand_per_capita = demand_per_capita(1:8784, 1:2);
% I hate NaN values
valid_demand_per_hour = demand_per_capita{:, 2};
valid_demand_per_hour(isnan(valid_demand_per_hour)) = 0;


 
total_demand_per_hour = valid_demand_per_hour * population;

% Load wind turbine power curve data
% Power in kW is stored in the second column
wind_power_curve = xlsread('turbine_power_curve_5_MW.xlsx', 'Sheet1', 'B2:E32');
wind_power_curve = wind_power_curve(:, 1:2);

% Load solar and wind data for each hour
solar_wind_data = readtable('solar_and_wind_data_hxh.csv', 'PreserveVariableNames', true);
original_solar_wind_data = solar_wind_data;

solar_wind_data = double(table2array(solar_wind_data(1:8784, 1:7)));

%% SOLAR PANEL CALCULATIONS
% Define efficiency and area parameters
solar_efficiency = 0.197 * 0.96; % cell + inverter
total_panel_area = 2.80; % m²
solar_irradiance = solar_wind_data(:, 5) / 1000; % Convert from Wh to kWh  x m2 

% Remove irradiance affected by shadows
solar_irradiance(solar_wind_data(:, 6) > 80) = 0;

% Compute energy output of a single panel
one_solar_energy_per_hour = solar_irradiance .* total_panel_area .* solar_efficiency;
total_energy_one_panel = sum(one_solar_energy_per_hour);

%% WIND TURBINE CALCULATIONS
% Adjust wind data for height and roughness
wind_data = solar_wind_data(:, 7);
h1 = 50; h2 = 140; z0 = 1.6;
wind_data = wind_data .* (log(h2 / z0) / log(h1 / z0));
wind_data = max(min(wind_data, max(wind_power_curve(:, 1))), min(wind_power_curve(:, 1))); % Interpolation limits

% Interpolate wind turbine power curve
one_wind_energy_per_hour = interp1(wind_power_curve(:, 1), wind_power_curve(:, 2), wind_data, 'linear', 'extrap');
total_energy_one_turbine = sum(one_wind_energy_per_hour);

%% HYDRO AND NUCLEAR POWER
% Define hydro and nuclear power outputs
hydropower = 20 * area; % in kW
hydropower_per_hour = repmat(hydropower, 8784, 1);

nuclear_power = 50000; % in kW
one_nuclear_power_per_hour = repmat(nuclear_power, 8784, 1);

%% OPTIMIZATION

best_chp_usage = inf;               %inizializzazione del chp usage vector
best_solution = [0, 0, 0, 0];       %inizializzazione dell' optimum vector

for n_pv = 355225:1:355225    %% after 5-6 calculation , I shortened the three ranges to speed up the final optimum computation
    for n_wind = 75:1:100
        for n_nuclear = 23:1:24
            for storage_capacity = 1207650:1:1207650

                % Compute CHP usage
                [total_chp_used, chp_usage_hourly, nuclear_usage_hourly, hydro_usage_hourly] = compute_chp_usage(n_pv, n_wind, n_nuclear, storage_capacity, ...
                    one_solar_energy_per_hour, one_wind_energy_per_hour, one_nuclear_power_per_hour, hydropower_per_hour, total_demand_per_hour);

                % Compute total cost
                total_cost = compute_cost(n_pv, n_wind, n_nuclear, storage_capacity);

                % Check if solution satisfies constraints and is optimal
                if total_cost <= max_investment && total_chp_used < best_chp_usage     %% if those constraints are validated, a new optimum replaces the previous
                    best_chp_usage = total_chp_used;
                    best_solution = [n_pv, n_wind, n_nuclear, storage_capacity];
                    best_chp_usage_hourly = chp_usage_hourly;
                end
            end
        end
    end
end

%% OPTIMIZED RESULTS
% Extract the best solution
n_pv = best_solution(1);
n_wind = best_solution(2);
n_nuclear = best_solution(3);
storage_capacity = best_solution(4);

% Calculate costs from function
[total_cost, cost_pv, cost_wind, cost_nuclear, cost_storage] = compute_cost(n_pv, n_wind, n_nuclear, storage_capacity);

%% SAVING RESULTS AND PLOTTING

% Save optimization results to a file
current_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
fileID = fopen('optimization_results.txt', 'a');

fprintf(fileID, '\n----------------------------------------\n');
fprintf(fileID, 'Results optimization logged on: %s\n', current_time);

fprintf(fileID, '\nOptimized Results:\n');
fprintf(fileID, 'Number of PV: %.2f\n', n_pv);
fprintf(fileID, 'Number of Wind Turbines: %.2f\n', n_wind);
fprintf(fileID, 'Number of Nuclear Plants: %.2f\n', n_nuclear);
fprintf(fileID, 'Storage Capacity: %.2f kWh\n', storage_capacity);

fprintf(fileID, '\nCosts:\n');
fprintf(fileID, 'Cost of PV: $%.2f\n', cost_pv);
fprintf(fileID, 'Cost of Wind: $%.2f\n', cost_wind);
fprintf(fileID, 'Cost of Nuclear: $%.2f\n', cost_nuclear);
fprintf(fileID, 'Cost of Storage: $%.2f\n', cost_storage);
fprintf(fileID, 'Total Cost: $%.2f\n', cost_pv + cost_wind + cost_nuclear + cost_storage);

fprintf(fileID, 'Total CHP used: %.2f kWh\n', best_chp_usage);
fclose(fileID);


% Plot CHP usage
figure;
plot(1:8784, best_chp_usage_hourly);
title('Hourly CHP Usage');
xlabel('Hour');
ylabel('CHP Usage (kWh)');
grid on;

% Plot energy production from PV and Wind
figure;
plot(1:8784, n_pv * one_solar_energy_per_hour, 1:8784, n_wind * one_wind_energy_per_hour);
title('Hourly Energy Production from PV and Wind');
legend('PV Energy', 'Wind Energy');
xlabel('Hour');
ylabel('Energy (kWh)');
grid on;

% Plot energy production from Nuclear and Hydro
figure;
plot(1:8784, n_nuclear*one_nuclear_power_per_hour, 1:8784, hydropower_per_hour);
title('Hourly Energy Production from Nuclear and Hydro');
legend('Nuclear Energy', 'Hydro Energy');
xlabel('Hour');
ylabel('Energy (kWh)');
grid on;

% Plot total energy demand and total production
total_energy_production = n_pv * one_solar_energy_per_hour + n_wind * one_wind_energy_per_hour + n_nuclear*one_nuclear_power_per_hour + hydropower_per_hour;
figure;
plot(1:8784, total_demand_per_hour, 1:8784, total_energy_production);
title('Hourly Total Energy Demand and Production');
legend('Total Demand', 'Total Production');
xlabel('Hour');
ylabel('Energy (kWh)');
grid on;

% Combined plot for CHP, PV, Wind, Nuclear, and Hydro usage as stacked layers using stacked bar chart
stacked_data = [n_nuclear*one_nuclear_power_per_hour,hydropower_per_hour, n_pv * one_solar_energy_per_hour, n_wind * one_wind_energy_per_hour,best_chp_usage_hourly,];

figure;
bar(1:8784, stacked_data, 'stacked', 'EdgeColor', 'none');

% Adding labels and legend
title('Hourly Energy Usage and Production (Stacked)');
legend('Nuclear Energy','Hydro Energy' ,'PV Energy', 'Wind Energy', 'CHP Usage');
xlabel('Hour');
ylabel('Energy (kWh)');
grid on;

% Plot total energy demand and total production
figure;
plot(1:8784, total_demand_per_hour, 'DisplayName', 'Total Demand');
hold on;
plot(1:8784, total_energy_production, 'DisplayName', 'Total Production');
hold off;

title('Hourly Total Energy Demand and Production');
legend;
xlabel('Hour');
ylabel('Energy (kWh)');
grid on;


%% FUNCTION DEFINITIONS
% Function to compute CHP usage
function [total_chp_used, chp_usage_hourly, nuclear_usage_hourly, hydro_usage_hourly] = compute_chp_usage(n_pv, n_wind, n_nuclear, storage_capacity, solar_energy_per_hour, wind_energy_per_hour, nuclear_power_per_hour, hydropower_per_hour, total_demand)
    energy_pv = n_pv * solar_energy_per_hour;
    energy_wind = n_wind * wind_energy_per_hour;
    energy_nuclear = n_nuclear * nuclear_power_per_hour;
    energy_hydro = hydropower_per_hour;

    energy_storage = zeros(8784, 1);
    energy_storage(1) = storage_capacity * 0.5; % Initial storage at 50%

    chp_usage_hourly = zeros(8784, 1);
    nuclear_usage_hourly = zeros(8784, 1);
    hydro_usage_hourly = zeros(8784, 1);

    for t = 1:8784
        energy_available = energy_pv(t) + energy_wind(t) + energy_nuclear(t) + energy_hydro(t);
        nuclear_usage_hourly(t) = energy_nuclear(t);  % (hydro and nuclear) constant
        hydro_usage_hourly(t) = energy_hydro(t);
        if energy_available >= total_demand(t)
            excess_energy = energy_available - total_demand(t);
            energy_storage(t + 1) = min(energy_storage(t) + 0.85 * excess_energy, storage_capacity);   % 0.85 è l'energy storage efficiency
        else
            residual_demand = total_demand(t) - energy_available;
            max_storage_use = min(energy_storage(t), storage_capacity / 6);   % decides between capacity remained in the battery and max speed 
            if residual_demand <= max_storage_use
                energy_storage(t + 1) = energy_storage(t) - residual_demand;
            else
                energy_storage(t + 1) = energy_storage(t) - max_storage_use;
                residual_demand = residual_demand - max_storage_use;
                chp_usage_hourly(t) = residual_demand;
            end
        end
    end
    total_chp_used = sum(chp_usage_hourly(~isnan(chp_usage_hourly)));     %no NAN values
end

% Function to compute total cost
function [total_cost, cost_pv, cost_wind, cost_nuclear, cost_storage] = compute_cost(n_pv, n_wind, n_nuclear, storage_capacity)
    % Calculate individual costs
    cost_pv = n_pv * 600;
    cost_wind = n_wind * (6.5e6);
    cost_nuclear = n_nuclear * (250e6);
    cost_storage = storage_capacity * 100;

    % Calculate total cost
    total_cost = cost_pv + cost_wind + cost_nuclear + cost_storage;
end




