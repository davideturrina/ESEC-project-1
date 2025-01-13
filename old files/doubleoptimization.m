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


% I hate NaN values
valid_demand_per_hour = demand_per_capita{:, 2};
valid_demand_per_hour(isnan(valid_demand_per_hour)) = 0;


total_demand_per_hour = valid_demand_per_hour * population;







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
fprintf('Total energy for one PV (2.80 m²) in a year: %.2f kWh\n', total_energy_one_panel);

%% Energy output of a single turbine
wind_speeds = wind_power_curve(:, 1);
power_output = wind_power_curve(:, 2);
wind_energy_per_hour = interp1(wind_speeds, power_output, wind_data, 'linear', 'extrap');
total_energy_one_turbine = sum(wind_energy_per_hour);
fprintf('Total energy produced for one turbine: %.2f kWh\n', total_energy_one_turbine);

%% Hydro and nuclear power
hydropower = 20 * area; % in kW
hydropower_per_hour = repmat(hydropower, 8784, 1);

nuclear_power = 50000; % in kW
nuclear_power_per_hour = repmat(nuclear_power, 8784, 1);

%% Initial storage setup
energy_storage = zeros(8785, 1);
energy_storage(1) = storage_capacity * 0.5; % Initial storage at 50%.

%% START OPTIMIZATION
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% Define optimization variables
x0 = [n_pv, n_wind, n_nuclear, storage_capacity]; % Initial guesses

% Bounds for variables
lb = [0, 0, 0, 0]; % No negative components
ub = [inf, inf , inf , inf]; % Arbitrary high bounds

% Objective function (minimize CHP usage)
objective = @(x) compute_chp_usage(x(1), x(2), x(3), x(4), solar_energy_per_hour, wind_energy_per_hour, nuclear_power_per_hour, hydropower_per_hour, total_demand_per_hour);

% Constraint function (cost <= max_investment)
constraints = @(x) deal([], compute_cost(x(1), x(2), x(3), x(4)) - max_investment);

% Optimization
[x_opt, fval] = fmincon(objective, x0, [], [], [], [], lb, ub, constraints, options);

% Optimized results
n_pv = ceil(x_opt(1));
n_wind = ceil(x_opt(2));
n_nuclear = floor(x_opt(3));                %questo è approssimato per difetto perchè altrimenti sforerebbe il budjet
storage_capacity = ceil(x_opt(4));


% Costs
cost_pv = n_pv * 600;
cost_wind = n_wind * (6.5e6);
cost_nuclear = n_nuclear * (250e6);
cost_storage = storage_capacity * 100;
cost_total= cost_pv+ cost_wind+ cost_nuclear+ cost_storage;


%% Optimized energy storage profile
energy_pv_opt = n_pv * solar_energy_per_hour;
energy_wind_opt = n_wind * wind_energy_per_hour;
energy_nuclear_opt = n_nuclear * nuclear_power_per_hour;
energy_hydro_opt = hydropower_per_hour;

energy_storage_opt = zeros(8784, 1);
energy_storage_opt(1) = storage_capacity * 0.5; % Initial storage at 50%

chp_usage_opt = zeros(8784, 1);

for t = 1:8784
    energy_available = energy_pv_opt(t) + energy_wind_opt(t) + energy_nuclear_opt(t) + energy_hydro_opt(t);
    if energy_available >= total_demand_per_hour(t)
        excess_energy = energy_available - total_demand_per_hour(t);
        energy_storage_opt(t + 1) = min(energy_storage_opt(t) + 0.85 * excess_energy, storage_capacity);
    else
        residual_demand = total_demand_per_hour(t) - energy_available;
        max_storage_use = min(energy_storage_opt(t), storage_capacity / 6);
        if residual_demand <= max_storage_use
            energy_storage_opt(t + 1) = energy_storage_opt(t) - residual_demand;
        else
            energy_storage_opt(t + 1) = energy_storage_opt(t) - max_storage_use;
            residual_demand = residual_demand - max_storage_use;
            chp_usage_opt(t) = residual_demand;
        end
    end
end

% Plot optimized energy storage
figure;
plot(energy_storage_opt, 'DisplayName', 'Optimized Storage', 'LineWidth', 1.5, 'Color', '#0072BD');
grid on;
xlabel('Hours');
ylabel('Stored Energy (kWh)');
legend('Location', 'best');
title('Energy Storage Profile');

% Energy contributions plot
figure;
plot(energy_pv_opt, 'DisplayName', 'PV', 'LineWidth', 1.5, 'Color', '#FF5733');
hold on;
plot(energy_wind_opt, 'DisplayName', 'Wind', 'LineWidth', 1.5, 'Color', '#33FF57');
plot(energy_hydro_opt, 'DisplayName', 'Hydro', 'LineWidth', 1.5, 'Color', '#3357FF');
hold off;
grid on;
xlabel('Hours');
ylabel('Energy (kWh)');
legend('Location', 'best');
title('Renewable Energy');


% Ottieni l'ora e la data correnti
current_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');

% Apri il file in modalità append
fileID = fopen('optimization_results.txt', 'a');

% Scrivi i risultati nel file con un'intestazione temporale
fprintf(fileID, '\n----------------------------------------\n');
fprintf(fileID, 'Results untitled4 logged on: %s\n', current_time);
fprintf(fileID, 'Optimized Results:\n');
fprintf(fileID, '- Number of PV panels: %d\n', n_pv);
fprintf(fileID, '- Number of wind turbines: %d\n', n_wind);
fprintf(fileID, '- Number of nuclear plants: %d\n', n_nuclear);
fprintf(fileID, '- Storage capacity: %d kWh\n', storage_capacity);

fprintf(fileID, '\nCosts:\n');
fprintf(fileID, '- Cost of PV: $%.2f\n', cost_pv);
fprintf(fileID, '- Cost of Wind: $%.2f\n', cost_wind);
fprintf(fileID, '- Cost of Nuclear: $%.2f\n', cost_nuclear);
fprintf(fileID, '- Cost of Storage: $%.2f\n', cost_storage);
fprintf(fileID, '- Total Cost: $%.2f\n', cost_pv + cost_wind + cost_nuclear + cost_storage);

fclose(fileID);

% CHP Usage Plot
figure;
plot(chp_usage_opt, 'DisplayName', 'CHP Usage', 'LineWidth', 1.5, 'Color', '#D95319');
grid on;
xlabel('Hours');
ylabel('CHP Usage (kWh)');
legend('Location', 'best');
title('CHP Usage Profile');

% Controlla se ci sono NaN in chp_usage_opt
if any(isnan(chp_usage_opt))
    warning('chp_usage_opt contiene valori NaN. Verranno ignorati nel calcolo della somma.');
    chp_usage_opt = chp_usage_opt(~isnan(chp_usage_opt));
end

% Total CHP Usage
total_chp_usage = sum(chp_usage_opt);
fprintf('Total CHP usage in the year: %.2f kWh\n', total_chp_usage);

% Log CHP Usage in File
fileID = fopen('optimization_results.txt', 'a');
fprintf(fileID, '\nCHP Usage:\n');
fprintf(fileID, '- Total CHP Usage: %.2f kWh\n', total_chp_usage);
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
