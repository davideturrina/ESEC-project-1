close all; 
clear;
clc;

%% Variables to find
capacity_storage = 0;
n_pv = 0;        
n_wind = 0;      
n_nuclear = 0;


%% DATA EXTRACTION
population = 2273800;
Area = 11246.8; % km²
MaxInvestment = population * 3000;

demandxcapita = readtable('demand_data_hxh_8784h.csv', 'PreserveVariableNames', true);
demandxcapita = demandxcapita(1:8784, 1:2);

% I hate NaN values
valid_demand_per_hour = demandxcapita{:, 2};
valid_demand_per_hour(isnan(valid_demand_per_hour)) = 0;

demand_total_xhour = valid_demand_per_hour * population;

% Wind turbine power curve data (Power in kW: second column)
wind_power_curve = table2array(readtable('turbine_power_curve_5_MW.xlsx', 'Sheet', 'Sheet1', 'Range', 'B2:F32'));
wind_power_curve = wind_power_curve(:, 1:5);

% Solar and wind data for each hour in Orlando
solar_wind_data = readtable('solar_and_wind_data_hxh.csv', 'PreserveVariableNames', true);
solar_wind_data = double(table2array(solar_wind_data(1:8784, 1:7)));

% Solar efficiency and panel area
solar_efficiency = 0.197 * 0.96; % cell+inverter
panel_area = 2.80; % m²
solar_data5 = solar_wind_data(:, 5) / 1000; % Solar irradiance was in Wh, now is in kWh

solar_data5(solar_wind_data(:,6) > 80) = 0; % Accounting for shadows beyond 90-10 degrees

% Extract numerical values for solar and wind data
wind_data_input = solar_wind_data(:, 7);
% Accounting for roughness
h1 = 50; h2 = 140; z0 = 1.6;
wind_data_input = wind_data_input .* (log(h2 / z0) / log(h1 / z0)); 

%% Energy output single panel
one_solar_energy_per_hour = solar_data5 .* panel_area .* solar_efficiency;
one_panel_total_energy = sum(one_solar_energy_per_hour);
disp(['Total energy for one PV (2.80 m²) in a year: ', num2str(one_panel_total_energy), ' kWh']);

%% Energy output single turbine (with interpolation)
wind_speeds = wind_power_curve(:, 1); 
power_output = wind_power_curve(:, 5); 
one_wind_energy_per_hour = interp1(wind_speeds, power_output, wind_data_input, 'linear', 'extrap');
one_turbine_total_energy = sum(one_wind_energy_per_hour);
disp(['Total energy produced for one turbine: ', num2str(one_turbine_total_energy), ' kWh']);

%% Hydro and nuclear
Hydropower = 20 * Area; % in kW
Hydropower_per_hour = repmat(Hydropower, 8784, 1);

one_Nuclear = 50000; % in kW
one_nuclear_per_hour = repmat(one_Nuclear, 8784, 1);


%% START OPTIMIZATION
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% Define optimization variables
x0 = [n_pv, n_wind, n_nuclear, capacity_storage]; % Initial guesses

% Lower and upper bounds
lb = [0, 0, 0, 0]; % No negative components
ub = [inf, inf, inf, inf]; % Arbitrary high bounds

% Objective function (minimize CHP usage)
objective = @(x) compute_chp_usage(x(1), x(2), x(3), x(4), one_solar_energy_per_hour, one_wind_energy_per_hour, one_nuclear_per_hour, Hydropower_per_hour, demand_total_xhour);

% Constraint function (cost <= MaxInvestment)
constraints = @(x) deal([], compute_cost(x(1), x(2), x(3), x(4)) - MaxInvestment);

% Optimization
[x_opt, fval] = fmincon(@(x) objective(x), x0, [], [], [], [], lb, ub, constraints, options);

% Compute results with optimized values
[n_pv, n_wind, n_nuclear, capacity_storage] = deal(x_opt(1), x_opt(2), x_opt(3), x_opt(4));

[total_chp_used, energy_storage] = compute_chp_usage(n_pv, n_wind, n_nuclear, capacity_storage, one_solar_energy_per_hour, one_wind_energy_per_hour, one_nuclear_per_hour, Hydropower_per_hour, demand_total_xhour);

% Save energy storage to workspace
assignin('base', 'energy_storage', energy_storage);

% Display results
fprintf('Optimized values:\n');
fprintf('Number of PV panels: %d\n', n_pv);
fprintf('Number of wind turbines: %d\n', n_wind);
fprintf('Number of nuclear plants: %d\n', n_nuclear);
fprintf('Storage capacity: %d kWh\n', capacity_storage);

% Display costs
cost_pv = n_pv * 600;
cost_wind = n_wind * (6.5e6);
cost_nuclear = n_nuclear * (250e6);
cost_storage = capacity_storage * 100;

fprintf('Cost of PV: $%.2f\n', cost_pv);
fprintf('Cost of Wind: $%.2f\n', cost_wind);
fprintf('Cost of Nuclear: $%.2f\n', cost_nuclear);
fprintf('Cost of Storage: $%.2f\n', cost_storage);
fprintf('Total Cost: $%.2f\n', cost_pv + cost_wind + cost_nuclear + cost_storage);

fprintf('Total CHP used: %.2f kWh\n', fval);

%% FUNCTIONS

function [total_cph_used, energy_storage] = compute_chp_usage(n_pv, n_wind, n_nuclear, capacity_storage, one_solar_energy_per_hour, one_wind_energy_per_hour, nuclear_per_hour, Hydropower_per_hour, demand_total)
    energy_pv = n_pv * one_solar_energy_per_hour; 
    energy_wind = n_wind * one_wind_energy_per_hour;
    energy_nuclear = n_nuclear * nuclear_per_hour;
    energy_hydro = Hydropower_per_hour;

    energy_storage = zeros(8784, 1);
    energy_storage(1) = 1.20269995e+06 * 0.5; % Initial storage at 50%.

    chp_usage = zeros(8784, 1);
    
    for t = 1:8784
        energy_available = energy_pv(t) + energy_wind(t) + energy_nuclear(t) + energy_hydro(t);
        
        if energy_available >= demand_total(t)
            excess_energy = energy_available - demand_total(t);
            energy_storage(t + 1) = min(energy_storage(t) + 0.85 * excess_energy, capacity_storage);
        else
            residual_demand = demand_total(t) - energy_available;
            max_storage_use = min(energy_storage(t), capacity_storage / 6);
            
            if residual_demand <= max_storage_use
                energy_storage(t + 1) = energy_storage(t) - residual_demand;
            else
                energy_storage(t + 1) = energy_storage(t) - max_storage_use;
                residual_demand = residual_demand - max_storage_use;
                chp_usage(t) = residual_demand;
            end
        end
    end

    total_cph_used = sum(chp_usage(~isnan(chp_usage)));
    % Truncate extra value in energy_storage to match demand
    energy_storage = energy_storage(1:8784);
end

function total_cost = compute_cost(n_pv, n_wind, n_nuclear, capacity_storage)
    cost_pv = n_pv * 600;
    cost_wind = n_wind * (6.5e6);
    cost_nuclear = n_nuclear * (250e6);
    cost_storage = capacity_storage * 100;
    total_cost = cost_pv + cost_wind + cost_nuclear + cost_storage;
end
