close all;
clear;
clc;

%% Variables to find
storage_capacity = 0;
n_pv = 0;
n_wind = 0;
n_nuclear = 0;

%% DATA EXTRACTION
population = 2273800;
area = 11246.8; % km²
max_investment = population * 3000;

demand_per_capita = readtable('demand_data_hxh_8784h.csv', 'PreserveVariableNames', true);
demand_per_capita = demand_per_capita(1:8784, 1:2);
total_demand_per_hour = demand_per_capita{:, 2} * population;

% Wind turbine power curve data (Power in kW: second column)
wind_power_curve = xlsread('turbine_power_curve_5_MW.xlsx', 'Sheet1', 'B2:E32');
wind_power_curve = wind_power_curve(:, 1:2);

% Solar and wind data for each hour in Orlando
solar_wind_data = readtable('solar_and_wind_data_hxh.csv', 'PreserveVariableNames', true);
solar_wind_data = double(table2array(solar_wind_data(1:8784, 1:7)));

% Solar efficiency and panel area
solar_efficiency = 0.197 * 0.96; % cell+inverter
panel_area = 2.80; % m²
solar_irradiance = solar_wind_data(:, 5) / 1000; % Convert from Wh to kWh

% Remove irradiance affected by shadows
solar_irradiance(solar_wind_data(:, 6) > 80) = 0;

% Adjust wind data for height and roughness
wind_data = solar_wind_data(:, 7);
h1 = 50; h2 = 140; z0 = 1.6;
wind_data = wind_data .* (log(h2 / z0) / log(h1 / z0));
wind_data = max(min(wind_data, max(wind_power_curve(:, 1))), min(wind_power_curve(:, 1))); % Interpolation limits

%% Energy output of a single panel
solar_energy_per_hour = solar_irradiance .* panel_area .* solar_efficiency;
total_energy_one_panel = sum(solar_energy_per_hour);
%%fprintf('Total energy for one PV (2.80 m²) in a year: %.2f kWh\n', total_energy_one_panel);

%% Energy output of a single turbine
wind_speeds = wind_power_curve(:, 1);
power_output = wind_power_curve(:, 2);
wind_energy_per_hour = interp1(wind_speeds, power_output, wind_data, 'linear', 'extrap');

total_energy_one_turbine = sum(wind_energy_per_hour);
%%fprintf('Total energy produced for one turbine: %.2f kWh\n', total_energy_one_turbine);

%% Hydro and nuclear power
hydropower = 20 * area; % in kW
hydropower_per_hour = repmat(hydropower, 8784, 1);

nuclear_power = 50000; % in kW
nuclear_power_per_hour = repmat(nuclear_power, 8784, 1);

%% Initial heuristic optimization
best_chp_usage = inf;
best_solution = [0, 0, 0, 0];

for n_pv = 402200:1:402300       % optimum [402229, 71, 24, 1185626] chp usage 818038379.82 kwh
    for n_wind = 70:1:71
        for n_nuclear =23:1:24
            for storage_capacity = 1185500:1:1185720
                % Compute CHP usage
                total_chp_used = compute_chp_usage(n_pv, n_wind, n_nuclear, storage_capacity, ...
                    solar_energy_per_hour, wind_energy_per_hour, nuclear_power_per_hour, hydropower_per_hour, total_demand_per_hour);

                % Compute total cost
                total_cost = compute_cost(n_pv, n_wind, n_nuclear, storage_capacity);

                % Check investment constraint
                if total_cost <= max_investment && total_chp_used < best_chp_usage
                    best_chp_usage = total_chp_used;
                    best_solution = [n_pv, n_wind, n_nuclear, storage_capacity];
                end
            end
        end
    end
end

% Optimized results
n_pv = best_solution(1);
n_wind = best_solution(2);
n_nuclear = best_solution(3);
storage_capacity = best_solution(4);

% Costs
cost_pv = n_pv * 600;
cost_wind = n_wind * (6.5e6);
cost_nuclear = n_nuclear * (250e6);
cost_storage = storage_capacity * 100;

% Get the current date and time
current_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');

% Open the file in append mode
fileID = fopen('optimization_results.txt', 'a');

% Write the results to the file with a timestamp header
fprintf(fileID, '\n----------------------------------------\n');
fprintf(fileID, 'Results untitled5 logged on: %s\n', current_time);

% Write optimized results
fprintf(fileID, '\nOptimized Results:\n');
fprintf(fileID, 'Number of PV: %.2f\n', n_pv);
fprintf(fileID, 'Number of Wind Turbines: %.2f\n', n_wind);
fprintf(fileID, 'Number of Nuclear Plants: %.2f\n', n_nuclear);
fprintf(fileID, 'Storage Capacity: %.2f kWh\n', storage_capacity);

% Write costs
fprintf(fileID, '\nCosts:\n');
fprintf(fileID, 'Cost of PV: $%.2f\n', cost_pv);
fprintf(fileID, 'Cost of Wind: $%.2f\n', cost_wind);
fprintf(fileID, 'Cost of Nuclear: $%.2f\n', cost_nuclear);
fprintf(fileID, 'Cost of Storage: $%.2f\n', cost_storage);
fprintf(fileID, 'Total Cost: $%.2f\n', cost_pv + cost_wind + cost_nuclear + cost_storage);

% Write total CHP usage
fprintf(fileID, 'Total CHP used: %.2f kWh\n', best_chp_usage);

% Close the file
fclose(fileID);


%% Function definitions
function total_chp_used = compute_chp_usage(n_pv, n_wind, n_nuclear, storage_capacity, solar_energy_per_hour, wind_energy_per_hour, nuclear_power_per_hour, hydropower_per_hour, total_demand)
    energy_pv = n_pv * solar_energy_per_hour;
    energy_wind = n_wind * wind_energy_per_hour;
    energy_nuclear = n_nuclear * nuclear_power_per_hour;
    energy_hydro = hydropower_per_hour;

    energy_storage = zeros(8784, 1);
    energy_storage(1) = storage_capacity * 0.5; % Initial storage at 50%

    chp_usage = zeros(8784, 1);
    for t = 1:8784
        energy_available = energy_pv(t) + energy_wind(t) + energy_nuclear(t) + energy_hydro(t);
        if energy_available >= total_demand(t)
            excess_energy = energy_available - total_demand(t);
            energy_storage(t + 1) = min(energy_storage(t) + 0.85 * excess_energy, storage_capacity);
        else
            residual_demand = total_demand(t) - energy_available;
            max_storage_use = min(energy_storage(t), storage_capacity / 6);
            if residual_demand <= max_storage_use
                energy_storage(t + 1) = energy_storage(t) - residual_demand;
            else
                energy_storage(t + 1) = energy_storage(t) - max_storage_use;
                residual_demand = residual_demand - max_storage_use;
                chp_usage(t) = residual_demand;
            end
        end
    end
    total_chp_used = sum(chp_usage(~isnan(chp_usage)));
end

function total_cost = compute_cost(n_pv, n_wind, n_nuclear, storage_capacity)
    cost_pv = n_pv * 600;
    cost_wind = n_wind * (6.5e6);
    cost_nuclear = n_nuclear * (250e6);
    cost_storage = storage_capacity * 100;

    total_cost = cost_pv + cost_wind + cost_nuclear + cost_storage;
end
