close all; 
clear;   
clc;
%% DATA EXTRACTION

                                                                                    % demand h - kWH 
demand = 'demand_data_hxh_8784h.csv';
demand = readtable(demand, 'PreserveVariableNames', true);
demand = demand(1:8784, 1:2);


wind_powercurve = xlsread('turbine_power_curve_5_MW.xlsx','Sheet1','B2:E32'); %Power in KW : second column

solar_wind_hxh= 'solar_and_wind_data_hxh.csv'; 
solar_wind_hxh= readtable(solar_wind_hxh, 'PreserveVariableNames', true);
solar_wind_hxh= solar_wind_hxh(1:8784, 1:7);                                  %solar 5row and wind 7row data in Orlando 
solar_efficiency = 0.197;
Area=2.80;



sol5_numeric = double(solar_wind_hxh{:,5});
wind7_numeric = double(solar_wind_hxh{:, 7});





%% VARIOUS PLOTS

% Plot Demand Data
figure;
plot(demand{:,1}, demand{:,2}, 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Demand (kW)');
title('Energy Demand');
grid on;

% Plot Wind Power Curve
figure;
plot(wind_powercurve(:,1), wind_powercurve(:,2), 'LineWidth', 1.5);
xlabel('Wind Speed (m/s)');
ylabel('Power Output (kW)');
title('Wind Turbine Power Curve');
grid on;

% Plot Solar and Wind Data con due assi y
time = 1:height(solar_wind_hxh); % Creazione dell'asse temporale (ore)
figure;

% Asse Y per i dati solari
yyaxis left;
plot(time, solar_wind_hxh{:,5}, 'DisplayName', 'Solar Data', 'LineWidth', 1.5);
ylabel('Solar Power (kW)');
ylim([0 max(solar_wind_hxh{:,5})*1.1]); % Margine superiore del 10%

% Asse Y per i dati eolici
yyaxis right;
plot(time, solar_wind_hxh{:,7}, 'DisplayName', 'Wind Data', 'LineWidth', 1.5);
ylabel('Wind Power (kW)');
ylim([0 max(solar_wind_hxh{:,7})*1.1]); % Margine superiore del 10%

% Titoli e asse X
xlabel('Time (hours)');
title('Solar and Wind Power Data');
legend('Solar Power', 'Wind Power');
grid on;

%% POWER ANALYSIS



OnePanelEnergyOUTPUT = sol5_numeric .* Area .* solar_efficiency;

figure; % Crea una nuova finestra per il grafico
plot(OnePanelEnergyOUTPUT);
xlabel('Index');
ylabel('Energy (E)');
title('Energy Output (E) Plot');
grid on;


%1turbineoutput...working on this



