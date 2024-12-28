close all;
clear;
clc;


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
wind_power_curve = table2array(readtable('turbine_power_curve_5_MW.xlsx', 'Sheet', 'Sheet1', 'Range', 'B2:F32'));

wind_power_curve = wind_power_curve(:, 1:5);

% Load solar and wind data for each hour
solar_wind_data = readtable('solar_and_wind_data_hxh.csv', 'PreserveVariableNames', true);
original_solar_wind_data = solar_wind_data;


solar_wind_data = double(table2array(solar_wind_data(1:8784, 1:7)));

%% SOLAR PANEL CALCULATIONS
% Define efficiency and area parameters
solar_efficiency = 0.197 * 0.96; % cell + inverter
total_panel_area = 2.80; % m²
solar_irradiance = solar_wind_data(:, 5) / 1000; % Convert from Wh to kWh       x m2 

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

%% Interpolate wind turbine power curve
one_wind_energy_per_hour = interp1(wind_power_curve(:, 1), wind_power_curve(:, 5), wind_data, 'linear', 'extrap');
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

for n_pv = 350730:1:350758    %% after 5-6 calculation , I shortened the three ranges to speed up the final optimum computation 24 and 75 are set to find pv and storage

    for n_wind = 75:1:75
        for n_nuclear = 24:1:24
            for storage_capacity = 1234540:1:1234586

                % Compute CHP usage
                [total_chp_used, chp_usage_hourly, storage_usage_hourly, energy_storage] = compute_chp_usage(n_pv, n_wind, n_nuclear, storage_capacity, ...
                    one_solar_energy_per_hour, one_wind_energy_per_hour, one_nuclear_power_per_hour, hydropower_per_hour, total_demand_per_hour);

                % Compute total cost
                total_cost = compute_cost(n_pv, n_wind, n_nuclear, storage_capacity);

                % Check if solution satisfies constraints and is optimal
                if total_cost <= max_investment && total_chp_used < best_chp_usage     %% if those constraints are validated, a new optimum replaces the previous
                    best_chp_usage = total_chp_used;
                    best_solution = [n_pv, n_wind, n_nuclear, storage_capacity];
                    best_chp_usage_hourly = chp_usage_hourly;
                    best_storage_usage_hourly=storage_usage_hourly;
                    best_energy_storage=energy_storage;
                    
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

%calculate storage used

total_storage_used=sum(best_storage_usage_hourly);

%calcalate energy mix


total_nuclear=sum(n_nuclear*one_nuclear_power_per_hour);
total_hydro=hydropower*8784;
total_solar=sum(n_pv * one_solar_energy_per_hour);
total_wind=sum(n_wind * one_wind_energy_per_hour);
total_chp_used = sum(chp_usage_hourly(~isnan(chp_usage_hourly)));     %no NAN values

total_energy_production = n_nuclear*one_nuclear_power_per_hour+hydropower_per_hour+ n_pv * one_solar_energy_per_hour+ n_wind * one_wind_energy_per_hour+best_chp_usage_hourly;



%% START OF PLOTS AND RESULTS 

%% STACKPILE

stacked_data = [n_nuclear*one_nuclear_power_per_hour,hydropower_per_hour, n_pv * one_solar_energy_per_hour, n_wind * one_wind_energy_per_hour,best_chp_usage_hourly,best_storage_usage_hourly];
% Define the start and end hours for the central week of each month
months_central_days = [15, 15+7]; % Days for the central week (middle of the month)
months_start_hours = (0:11) * 30 * 24 + (months_central_days(1) - 1) * 24 + 1; % Start hours for each central week
months_end_hours = (0:11) * 30 * 24 + (months_central_days(2) - 1) * 24;       % End hours for each central week

% Iterate over each month and create a stacked bar plot
for month = 1:length(months_start_hours)
    % Define start and end hour for the current central week
    start_hour = months_start_hours(month);
    end_hour = months_end_hours(month);

    % Ensure the end_hour does not exceed the data size
    if end_hour > size(stacked_data, 1)
        end_hour = size(stacked_data, 1);
    end

    % Extract data for the selected central week
    week_data = stacked_data(start_hour:end_hour, :);

    % Create the stacked bar plot
    figure(month); % Assign unique figure number
    bar(1:size(week_data, 1), week_data, 'stacked', 'EdgeColor', 'none');

    % Add title, labels, and legend
    title(['Energy Usage and Production (Central Week of Month ' num2str(month) ')']);
    legend('Nuclear Energy', 'Hydro Energy', 'PV Energy', 'Wind Energy', 'CHP Usage', 'Storage Usage', 'Location', 'south');
    xlabel('Hour of the Week');
    ylabel('Energy (kWh)');
    grid on;
end

%% PIE

energy_sources = [total_chp_used, total_nuclear, total_solar, total_wind,total_hydro];
source_labels = {'CHP', 'Nuclear', 'Solar', 'Wind', 'Hydro'};

total_energy = sum(energy_sources);
percentages = (energy_sources / total_energy) * 100;

labels_with_percentages = cell(size(source_labels));
for i = 1:length(source_labels)
    labels_with_percentages{i} = sprintf('%s (%.1f%%)', source_labels{i}, percentages(i));
end

figure(length(months_start_hours) + 1); % Assign unique figure number
pie(energy_sources, labels_with_percentages);
title('Energy Sources Distribution', 'FontSize', 14);

colormap(jet);
set(gca, 'FontSize', 12);

%% STORAGE CAPACITY

start_hour = 1912; % Third week of March 20-27
end_hour = 2080;

march_week_storage = best_energy_storage(start_hour:end_hour);

time = 1:(end_hour - start_hour + 1);
figure(length(months_start_hours) + 2); % Assign unique figure number
plot(time, march_week_storage, 'LineWidth', 1.5);
grid on;

xlabel('Time (Hours)');
ylabel('Storage Capacity');
title('Power Storage Capacity - Third Week of March');

%% USEFUL GRAPH

start_hour = 6577;
end_hour = 6744;

weekly_demand = total_demand_per_hour(start_hour:end_hour);
weekly_production = total_energy_production(start_hour:end_hour);
weekly_storage_usage = best_storage_usage_hourly(start_hour:end_hour);
week_storage_capacity = energy_storage(start_hour:end_hour);

figure(length(months_start_hours) + 3); % Assign unique figure number
plot(1:length(weekly_demand), weekly_demand, 'DisplayName', 'Total Demand', 'LineWidth', 2);
hold on;
plot(1:length(weekly_production), weekly_production, 'DisplayName', 'Total Production', 'LineWidth', 2); 
plot(1:length(weekly_storage_usage), weekly_storage_usage, 'DisplayName', 'Storage Usage', 'LineWidth', 2); 
plot(1:length(weekly_storage_usage), week_storage_capacity, 'DisplayName', 'Storage Capacity', 'LineWidth', 2); 

hold off;

title('First Week of October: Hourly Total Energy Demand and Production and Usage');
legend;
xlabel('Hour (First Week of October)');
ylabel('Energy (kWh)');
grid on;

%% LOAD FACTOR

% load factor = average load / max load

total_demand_per_hour;
mean_value = mean(best_chp_usage);
max_value = max(best_chp_usage);

load_factor = mean_value / max_value;

%% COSTS

costs = [total_cost, cost_pv, cost_wind, cost_nuclear, cost_storage];
sources = {'PV', 'Wind', 'Nuclear', 'Storage'};

shares = (costs(2:end) / costs(1)) * 100;

figure(length(months_start_hours) + 4); % Assign unique figure number
pie(shares, sources);

title('Cost Share by Source');

% From https://app.electricitymaps.com/map, Florida values in g/kWh  
co2nuke = 0.012; 
co2sol = 0.026;
co2chp = 0.560;
co2wind = 0.011;
co2hydro = 0.024;

emissions_hydro = total_hydro * co2hydro;
emissions_wind = total_wind * co2wind;
emissions_nuclear = total_nuclear * co2nuke;
emissions_solar = total_solar * co2sol;
emissions_chp = best_chp_usage * co2chp;
total_emissions = emissions_hydro + emissions_wind + emissions_nuclear + emissions_solar + emissions_chp;

%% PIE EMISSIONS

% Combine data into arrays
emissions_values = [emissions_hydro, emissions_wind, emissions_nuclear, emissions_solar, emissions_chp];
emissions_labels = {'Hydro', 'Wind', 'Nuclear', 'Solar', 'CHP'};

% Create the pie chart
figure(length(months_start_hours) + 5); % Assign unique figure number
pie(emissions_values, emissions_labels);
title('Emissions Contribution by Source (kg CO2)');

%% ELECTRICITY SPARED  (storage is not accounted in the total_energy_prod)

diff_electricity = total_energy_production - total_demand_per_hour;
electricity_spared = sum(diff_electricity(diff_electricity > 0));
percentage_spared = (electricity_spared / total_energy) * 100;



%% SAVING RESULTS 

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

fprintf(fileID, '\nEmissions:\n');
fprintf(fileID, 'Emissions from Hydro: %.2f kg CO2\n', emissions_hydro);
fprintf(fileID, 'Emissions from Wind: %.2f kg CO2\n', emissions_wind);
fprintf(fileID, 'Emissions from Nuclear: %.2f kg CO2\n', emissions_nuclear);
fprintf(fileID, 'Emissions from Solar: %.2f kg CO2\n', emissions_solar);
fprintf(fileID, 'Emissions from CHP: %.2f kg CO2\n', emissions_chp);
fprintf(fileID, 'Total Emissions: %.2f kg CO2\n', total_emissions);

fprintf(fileID, 'Total CHP used: %.2f kWh\n', best_chp_usage);

fprintf(fileID, 'load_factor for chp: %.2f\n', load_factor);

fprintf(fileID, 'electricity spared in kWh: %.2f\n', electricity_spared);
fprintf(fileID, 'percentage spared: %.2f\n', percentage_spared);


fclose(fileID);


%% FUNCTION DEFINITIONS

function [total_chp_used_notopt, chp_usage_hourly, storage_usage_hourly, energy_storage] = compute_chp_usage(n_pv, n_wind, n_nuclear, storage_capacity, solar_energy_per_hour, wind_energy_per_hour, nuclear_power_per_hour, hydropower_per_hour, total_demand)
    energy_pv = n_pv * solar_energy_per_hour;
    energy_wind = n_wind * wind_energy_per_hour;
    energy_nuclear = n_nuclear * nuclear_power_per_hour;
    energy_hydro = hydropower_per_hour;

    energy_storage = zeros(8784, 1);
    storage_usage_hourly = zeros(8784, 1);
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
                storage_usage_hourly(t) = residual_demand;
            else
                energy_storage(t + 1) = energy_storage(t) - max_storage_use;
                residual_demand = residual_demand - max_storage_use;
                storage_usage_hourly(t) = max_storage_use;
                chp_usage_hourly(t) = residual_demand;
            end
        end
    end
    total_chp_used_notopt = sum(chp_usage_hourly(~isnan(chp_usage_hourly)));     %no NAN values
    
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




