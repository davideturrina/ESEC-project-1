clear

% demand h - kWH 
demand = 'demand_data_hxh_8784h.csv';
demand = readtable(demand, 'PreserveVariableNames', true);
demand = demand(1:8784, 1:2);


wind_powercurve = xlsread('turbine_power_curve_5_MW.xlsx','Sheet1','B2:E32'); %Power in KW : second column







