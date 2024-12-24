%% Separate Micro Water Energy Nexus Optimization Model
% Jesus Silva Rodrigues | University of Houston Department of Electrical
% and Computer Engineering

% This model implements the same MWEN network, but optimizing each system
% individually, with no participation or interaction with each other.

clear
close all
clc
%% Datasets
% LOADING DATA
SOLAR_coast = readmatrix('SOLAR_GAL.CSV','Range','C5091:M5834');
SOLAR_city = readmatrix('SOLAR_HOU.CSV','Range','C5091:M5834');
WIND_coast = readmatrix('METEO_GAL.xlsx','Range','G5081:G5824');
WIND_city = readmatrix('METEO_HOU.xlsx','Range','G50333:G57772');
P_LOAD1 = 1000*readmatrix('P_LOAD.xlsx','Range','B4:Y7');%[kW]
W_LOAD1 = readmatrix('WATER.xlsx','Range','B4:Y7');%[gal/h]

% PROCESSING DATA
% Solar Irradiances
[GHI_coast,DNI_coast,DHI_coast] = SolarIR(SOLAR_coast);
[GHI_city,DNI_city,DHI_city] = SolarIR(SOLAR_city);

% Wind Speeds
WS_coast = [];
for c = 1:24:length(WIND_coast)
    WS_coast = [WS_coast WIND_coast(c:c+23)];
end
WS_coast_av = mean(WS_coast,2);

WS_city_m = [];
for c = 1:10:length(WIND_city)
    WS_city_m = [WS_city_m WIND_city(c:c+9)];
end
WS_city_m = mean(WS_city_m)';
WS_city = [];
for c = 1:24:length(WS_city_m)
    WS_city = [WS_city WS_city_m(c:c+23)];
end
WS_city_av = mean(WS_city,2);

% Grid Hourly Prices
load('GRID.mat');% 24-hr grid-prices vector

delta_t = 1;%time-step [h]

%% Setting Parameters
% The numeric values of each matrix can be changed if a different test case
% is desired

% GENERATORS
% Matrix dimensions: (num of MGs, num of gens on each MG) 
P_G_min1 = [100 150 0;100 150 260;250 0 0;0 0 0];%[kW]
P_G_max1 = [650 800 0;650 800 940;900 0 0;0 0 0];%[kW]
C_G1 = [0.33 0.28 0;0.33 0.28 0.23;0.26 0 0;0 0 0];%[$/kWh]
SU_G1 = [15.00 13.00 0;15.00 13.00 10.35;11.85 0 0;0 0 0];%[$]
NL_G1 = [11.00 8.70 0;11.00 8.70 7.40;8.45 0 0;0 0 0];%[$/h]
u_G_init1 = zeros(size(P_G_min1));

% ENERGY STORAGE
% Matrix dimensions: (num of MGs, num of ESS on each MG)
P_ES_limit1 = [300 300 300 300 300 300;725 725 725 725 0 0;300 300 300 300 300 0;725 725 725 725 725 0];%[kW]
EL_ES_min1 = [0 0 0 0 0 0;225 225 225 225 0 0;0 0 0 0 0 0;225 225 225 225 225 0];%[kWh]
EL_ES_max1 = [1660 1660 1660 1660 1660 1660;2500 2500 2500 2500 0 0;1660 1660 1660 1660 1660 0;2500 2500 2500 2500 2500 0];%[kWh]
EL_ES_init1 = EL_ES_min1;
eff_ES_c1 = [0.80 0.80 0.80 0.80 0.80 0.80;0.95 0.95 0.95 0.95 1 1;0.80 0.80 0.80 0.80 0.80 1;0.95 0.95 0.95 0.95 0.95 1];%[-]
eff_ES_d1 = [0.60 0.60 0.60 0.60 0.60 0.60;0.98 0.98 0.98 0.98 1 1;0.60 0.60 0.60 0.60 0.60 1;0.98 0.98 0.98 0.98 0.98 1];%[-]
CW_ES1 = 0.0707*ones(size(P_ES_limit1));%[gal/kWh]
CW_ES1(P_ES_limit1 == 0) = 0; % Energy storage water consumtion. Primarily for HESS
CW_ES1(2,:) = 0;% MWEN 2 and MWEN 4 have BESS
CW_ES1(4,:) = 0;

% NETWORK POWER TIE-LINE LIMIT
P_E_limit1 = 1400*ones(4,1);%[kW]

% WASTEWATER TREATMENT
%Reclaimed Wastewater
w_rec = 0.50;% Percentaje of water reclamed
% Wastewater treatment characteristics
W_WW_min1 = [20;50;20;50];%[gal/h]
W_WW_max1 = [720;1000;720;1000];%[gal/h]
WL_WW_max1 = [12000;18600;12000;18600];%[gal]
WL_WW_init1 = 0.04*WL_WW_max1;%[gal]
NL_WW1 = [75.00;75.00;75.00;75.00];%[$/h]
CW_WW1 = [365;382;365;382];%[gal/kWh]
eff_WW_pump = 0.84;%[-]

% CONTROLLABLE WATER TREATMENT
W_WT_min1 = [120;35;150;0];%[gal/h]
W_WT_max1 = [1500;525;1730;0];%[gal/h]
NL_WT1 = [17.20;25.5;15.00;0];%[$/h]
CW_WT1 = [96;4400;100;0];%[gal/kWh]
eff_WT_pump = 0.84;%[-]

% WATER STORAGE
W_ST_limit1 = [900;900;900;900];%[gal/h]
WL_ST_limit1 = [10209;10209;10209;10209];%[gal]
WL_ST_init1 = zeros(size(W_ST_limit1));%#[gal]
eff_ST_pump = 0.84;%[-]

% NETWORK WATER TIE-LINE LIMIT
W_E_limit1 = 980*ones(4,1);%[gal/h]

% SOLAR POWER
S_SP = [10*50;15*180;10*40;15*220];%[m^2]
eff_SP = [0.30;0.35;0.30;0.35];%[-]
tilt_a = 27*ones(size(S_SP));%[°]
collect_a = zeros(size(S_SP));%[°]
Lat_coast = readmatrix('SOLAR_GAL.csv','Range','E1:E1');%[°]
Lat_city = readmatrix('SOLAR_HOU.csv','Range','E1:E1');%[°]
n = 224;% A day in August
P_SP1 = zeros(size(P_LOAD1));
for m = 1:2:size(P_SP1,1)
    P_SP1(m,:) = SolarPow(DNI_coast,DHI_coast,GHI_coast,n,Lat_coast,S_SP(m),eff_SP(m),tilt_a(m),collect_a(m));
end
for m = 2:2:size(P_SP1,1)
    P_SP1(m,:) = SolarPow(DNI_city,DHI_city,GHI_city,n,Lat_city,S_SP(m),eff_SP(m),tilt_a(m),collect_a(m));
end

% WIND POWER
Prated_WT = [250 250 300 300;0 0 0 0;250 250 250 250;0 0 0 0];%[kW]
CutIn = [2.4 2.4 3 3;0 0 0 0;2.4 2.4 2.4 2.4;0 0 0 0;];%[m/s]
v_rated = [5.6 5.6 7.2 7.2;0 0 0 0;5.6 5.6 5.6 5.6;0 0 0 0];%[m/s]
CutOut = [10 10 13.6 13.6;0 0 0 0;10 10 10 10;0 0 0 0];%[m/s]
P_WP1 = zeros(size(P_LOAD1));%[kW]
for m = 1:2:size(Prated_WT,1)
    P_WP1(m,:) = WindPow(Prated_WT(m,:),WS_coast_av,CutIn(m,:),v_rated(m,:),CutOut(m,:));
end
for m = 2:2:size(Prated_WT,1)
    P_WP1(m,:) = WindPow(Prated_WT(m,:),WS_city_av,CutIn(m,:),v_rated(m,:),CutOut(m,:));
end

% GRID AND NETWORK EXCHANGE
% Power
C_grid_plus = Grid';%[$/kWh]
C_grid_minus = 0.8*C_grid_plus;%[$/kWh]
C_Np = 0.9*C_grid_plus;%[$/kWh]
Cp_shed = 1000;% Power load shed penalty
% Water
C_Wmain_plus = 0.006;%[$/gal]
C_Nw = 0.9*C_Wmain_plus;%[$/gal]
Cw_shed = 1000;% Water load shed penalty

% Allocating Space for final value variables
[M,T] = size(P_LOAD1);
G = size(P_G_max1,2);
ES = size(P_ES_limit1,2);

u_G_f = zeros(M,T,G);
v_G_f = zeros(M,T,G);
P_G_f = zeros(M,T,G);
P_grid_plus_f = zeros(M,T);
P_grid_minus_f = zeros(M,T);
P_ESc_f = zeros(M,T,ES);
P_ESd_f = zeros(M,T,ES);
EL_ES_f = zeros(M,T,ES);
W_ES_f = zeros(M,T,ES);
e_ESc_f = zeros(M,T,ES);
e_ESd_f = zeros(M,T,ES);
P_shed_f = zeros(M,T);
W_WW_f = zeros(M,T);
W_WT_f = zeros(M,T);
P_WT_f = zeros(M,T);
P_WW_f = zeros(M,T);
WL_WW_f = zeros(M,T);
WR_WW_f = zeros(M,T);
u_WT_f = zeros(M,T);
u_WW_f = zeros(M,T);
W_STc_f = zeros(M,T);
W_STd_f = zeros(M,T);
WL_ST_f = zeros(M,T);
sp_ST_f = zeros(M,T);
sv_ST_f = zeros(M,T);
W_main_plus_f = zeros(M,T);
P_WT_pump_f = zeros(M,T);
P_ST_pump_f = zeros(M,T);
P_WW_pump_f = zeros(M,T);
W_shed_f = zeros(M,T);
f_cost_f = zeros(M,1);


for s = 1:M

% Separating parameters for MWEN m
P_G_min = P_G_min1(s,:);
P_G_max = P_G_max1(s,:);
C_G = C_G1(s,:);
SU_G = SU_G1(s,:);
NL_G = NL_G1(s,:);
u_G_init = u_G_init1(s,:);
P_ES_limit = P_ES_limit1(s,:);
EL_ES_min = EL_ES_min1(s,:);
EL_ES_max = EL_ES_max1(s,:);
EL_ES_init = EL_ES_init1(s,:);
eff_ES_c = eff_ES_c1(s,:);
eff_ES_d = eff_ES_d1(s,:);
CW_ES = CW_ES1(s,:);
P_E_limit = P_E_limit1(s);
W_WW_min = W_WW_min1(s);
W_WW_max = W_WW_max1(s);
WL_WW_max = WL_WW_max1(s);
WL_WW_init = WL_WW_init1(s);
NL_WW = NL_WW1(s);
CW_WW = CW_WW1(s);
W_WT_min = W_WT_min1(s);
W_WT_max = W_WT_max1(s);
NL_WT = NL_WT1(s);
CW_WT = CW_WT1(s);
W_ST_limit = W_ST_limit1(s);
WL_ST_limit = WL_ST_limit1(s);
WL_ST_init = WL_ST_init1(s);
W_E_limit = W_E_limit1(s);
P_SP = P_SP1(s,:);
P_WP = P_WP1(s,:);
P_grid_limit = P_E_limit;
P_LOAD = P_LOAD1(s,:);
W_main_limit = W_E_limit;
W_LOAD = W_LOAD1(s,:);

fprintf("MicroWEN %d Optimization...",s)

%% Individual MG Optimization

cvx_begin
% VARIABLES
    % Power Variables
    variable u_G(G,T) binary
    variable v_G(G,T) binary
    variables P_G(G,T) P_grid_plus(1,T) P_grid_minus(1,T)
    variables P_ESc(ES,T) P_ESd(ES,T) EL_ES(ES,T) W_ES(ES,T)
    variable e_ESc(ES,T) binary
    variable e_ESd(ES,T) binary
    variable p_plus(1,T) binary
    variable p_minus(1,T) binary
    variable P_shed(1,T)
    % Water Variables
    variables W_WW(1,T) W_WT(1,T) P_WT(1,T) P_WW(1,T) WL_WW(1,T) WR_WW(1,T)
    variable u_WT(1,T) binary
    variable u_WW(1,T) binary
    variables W_STc(1,T) W_STd(1,T) WL_ST(1,T)
    variable sp_ST(1,T) binary
    variable sv_ST(1,T) binary
    variables W_main_plus(1,T) P_WT_pump(1,T) P_ST_pump(1,T) P_WW_pump(1,T)
    variable W_shed(1,T)
    
% OBJECTIVE FUNCTION
    % Power cost
    f_E = 0;
    for t = 1:1:T
        for g = 1:1:G
            f_E = f_E + SU_G(g)*v_G(g,t) + delta_t*(NL_G(g)*u_G(g,t)+C_G(g)*P_G(g,t)); 
        end
        f_E = f_E + C_grid_plus(t)*P_grid_plus(t) - C_grid_minus(t)*P_grid_minus(t) + Cp_shed*P_shed(t);
    end
    % Water cost
    f_W = 0;
    for t = 1:1:T
        f_W = f_W + delta_t*(NL_WW*u_WW(t)+NL_WT*u_WT(t)) + C_Wmain_plus*W_main_plus(t) + Cw_shed*W_shed(t);
    end
    minimize f_E + f_W
    
% CONSTRAINTS
    subject to
        %P_shed(1) == 0
        % Power Constraints
        for t = 1:1:T
            % Generators
            for g = 1:1:G
                P_G_min(g)*u_G(g,t) <= P_G(g,t) <= P_G_max(g)*u_G(g,t)
                if t == 1
                    v_G(g,t) >= u_G(g,t) - u_G_init(g)
                else
                    v_G(g,t) >= u_G(g,t) - u_G(g,t-1)
                end
            end
            % Energy Storage
            for b = 1:1:size(P_ES_limit,2)
                0 <= P_ESc(b,t) <= P_ES_limit(b)*e_ESc(b,t)
                0 <= P_ESd(b,t) <= P_ES_limit(b)*e_ESd(b,t)
                e_ESc(b,t) + e_ESd(b,t) <= 1
                if t == 1
                    EL_ES(b,t) == EL_ES_init(b) + delta_t*(eff_ES_c(b)*P_ESc(b,t) - P_ESd(b,t)/eff_ES_d(b))
                else
                    EL_ES(b,t) == EL_ES(b,t-1) + delta_t*(eff_ES_c(b)*P_ESc(b,t) - P_ESd(b,t)/eff_ES_d(b))
                end
                EL_ES_min(b) <= EL_ES(b,t) <= EL_ES_max(b)
                W_ES(b,t) == CW_ES(b)*P_ESc(b,t)
            end
            % MicroWEN Power Balance
            sum(P_G(:,t)) + sum(P_ESd(:,t)-P_ESc(:,t)) + P_grid_plus(t) - P_grid_minus(t) == P_LOAD(t) - P_SP(t) - P_WP(t) + P_WW(t) + P_WW_pump(t) + P_WT(t) + P_WT_pump(t) + P_ST_pump(t) - P_shed(t)
            % Main Grid
            0 <= P_grid_plus(t) <= P_grid_limit*p_plus(t)
            0 <= P_grid_minus(t) <= P_grid_limit*p_minus(t)
            p_plus(t) + p_minus(t) <= 1
        end
        % Water Constraints
        for t = 1:1:T
            % Wastewater Treatment
            W_WW_min*u_WW(t) <= W_WW(t) <= W_WW_max*u_WW(t)
            if t == 1
                WL_WW(t) == WL_WW_init + delta_t*(WR_WW(t)-W_WW(t))
            else
                WL_WW(t) == WL_WW(t-1) + delta_t*(WR_WW(t)-W_WW(t))
            end
            0 <= WL_WW(t) <= WL_WW_max
            if t == 1
                WR_WW(t) == w_rec*(W_LOAD(end))
            else
                WR_WW(t) == w_rec*(W_LOAD(t-1)-W_shed(t-1))
            end
            W_WW(t) == CW_WW*P_WW(t)
            % Regular Water Treatment
            W_WT_min*u_WT(t) <= W_WT(t) <= W_WT_max*u_WT(t)
            W_WT(t) == CW_WT*P_WT(t)
            P_WT(t) >= 0
            % Water Storage
            0 <= W_STc(t) <= W_ST_limit*sp_ST(t)
            0 <= W_STd(t) <= W_ST_limit*sv_ST(t)
            sp_ST(t) + sv_ST(t) <= 1
            if t == 1
                WL_ST(t) == WL_ST_init + delta_t*(W_STc(t)-W_STd(t))
            else
                WL_ST(t) == WL_ST(t-1) + delta_t*(W_STc(t)-W_STd(t))
            end
            0 <= WL_ST(t) <= WL_ST_limit
            % MicroWEN Water Balance
            W_WT(t) + W_STd(t) - W_STc(t) + W_main_plus(t) == W_LOAD(t) - W_WW(t) + sum(W_ES(:,t)) - W_shed(t)
            % Main Water System
            0 <= W_main_plus(t) <= W_main_limit
            % Wastewater Pump
            1000*eff_WW_pump*P_WW_pump(t) == 2.3177*W_WW(t)
            % Water Treatment Pump
            1000*eff_WT_pump*P_WT_pump(t) == 2.3177*W_WT(t)
            % Water Storage Pump
            1000*eff_ST_pump*P_ST_pump(t) == 2.3177*W_STc(t)
        end
        % Energy and Water Storage Reliability
        sum(EL_ES(:,T)) == 5*0.3*max(P_LOAD)
        WL_ST(T) == 5*0.3*max(W_LOAD)
        WL_WW(T) <= 0.04*WL_WW_max
        % Power and Water Load Shed
        P_shed >= 0
        W_shed >= 0
    
cvx_end

%% Objective Values

% Power cost
f_E = 0;
for t = 1:1:T
    for g = 1:1:G
        f_E = f_E + SU_G(g)*v_G(g,t) + delta_t*(NL_G(g)*u_G(g,t)+C_G(g)*P_G(g,t));
    end
    f_E = f_E + C_grid_plus(t)*P_grid_plus(t) - C_grid_minus(t)*P_grid_minus(t);
end

% Water cost
f_W = 0;
for t = 1:1:T
    f_W = f_W + delta_t*(NL_WW*u_WW(t)+NL_WT*u_WT(t)) + C_Wmain_plus*W_main_plus(t);
end

% Total cost
f_cost = f_E + f_W;

%% Retrieving m Values into Final Variables
for t = 1:T
    for g = 1:G
        u_G_f(s,t,g) = u_G(g,t);
        v_G_f(s,t,g) = v_G(g,t);
        P_G_f(s,t,g) = P_G(g,t);
    end
    P_grid_plus_f(s,t) = P_grid_plus(t);
    P_grid_minus_f(s,t) = P_grid_minus(t);
    for b = 1:ES
        P_ESc_f(s,t,b) = P_ESc(b,t);
        P_ESd_f(s,t,b) = P_ESd(b,t);
        EL_ES_f(s,t,b) = EL_ES(b,t);
        W_ES_f(s,t,b) = W_ES(b,t);
        e_ESc_f(s,t,b) = e_ESc(b,t);
        e_ESd_f(s,t,b) = e_ESd(b,t);
    end
    P_shed_f(s,t) = P_shed(t);
    W_WW_f(s,t) = W_WW(t);
    W_WT_f(s,t) = W_WT(t);
    P_WT_f(s,t) = P_WT(t);
    P_WW_f(s,t) = P_WW(t);
    WL_WW_f(s,t) = WL_WW(t);
    WR_WW_f(s,t) = WR_WW(t);
    u_WT_f(s,t) = u_WT(t);
    u_WW_f(s,t) = u_WW(t);
    W_STc_f(s,t) = W_STc(t);
    W_STd_f(s,t) = W_STd(t);
    WL_ST_f(s,t) = WL_ST(t);
    sp_ST_f(s,t) = sp_ST(t);
    sv_ST_f(s,t) = sv_ST(t);
    W_main_plus_f(s,t) = W_main_plus(t);
    P_WT_pump_f(s,t) = P_WT_pump(t);
    P_ST_pump_f(s,t) = P_ST_pump(t);
    P_WW_pump_f(s,t) = P_WW_pump(t);
    W_shed_f(s,t) = W_shed(t);
end
f_cost_f(s) = f_cost;% Individual operation costs

end
%% Output Results
fprintf("Optimal Operation Costs\n")
for m = 1:M
    fprintf("MWEN%d: %0.2f\n",m,f_cost_f(m))
end
fprintf("TOTAL: %0.2f.\n",sum(f_cost_f))

% Net Loads
figure(1)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;P_LOAD1(m,:)-P_SP1(m,:)-P_WP1(m,:)+P_WW_f(m,:)+P_WW_pump_f(m,:)+P_WT_f(m,:)+P_WT_pump_f(m,:)+P_ST_pump_f(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Power [kW]')
axis([0 T+1 min(min(P_LOAD1-P_WP1-P_SP1+P_WW_f+P_WW_pump_f+P_WT_f+P_WT_pump_f+P_ST_pump_f))-100 max(max(P_LOAD1-P_WP1-P_SP1+P_WW_f+P_WW_pump_f+P_WT_f+P_WT_pump_f+P_ST_pump_f))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northwest')
title('Power Net Demand')

% Energy Storage charging/discharging
figure(2)
Y = [];
Leg = [];
P_ESc_all = sum(P_ESc_f,3);
P_ESd_all = sum(P_ESd_f,3);
for m = 1:M
    Y = [Y;P_ESc_all(m,:)-P_ESd_all(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Power [kW]')
legend(Leg,'Location','Southwest')
axis([0 T+1 min(min(P_ESc_all-P_ESd_all))-100 max(max(P_ESc_all-P_ESd_all))+100])
xticks(0:2:T+1)
title('Energy Storage Input/Output')

% Energy Storage charge level
figure(3)
Y = [];
Leg = [];
EL_ES_all = sum(EL_ES_f,3);
for m = 1:M
    Y = [Y;EL_ES_all(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Charge Level [kWh]')
legend(Leg,'Location','Best')
axis([0 T+1 min(min(EL_ES_all))-1000 max(max(EL_ES_all))+1000])
xticks(0:2:T+1)
title('Energy Storage Charge Level')

% Generators Output
figure(4)
Y = [];
Leg = [];
P_G_all = sum(P_G_f,3);
for m = 1:M
    Y = [Y;P_G_all(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Power [kW]')
legend(Leg,'Location','Northwest')
axis([0 T+1 min(min(P_G_all))-100 max(max(P_G_all))+100])
xticks(0:2:T+1)
title('Generators')

% Main Grid
figure(5)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;P_grid_plus_f(m,:)-P_grid_minus_f(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
Y = [Y;sum(P_grid_plus_f-P_grid_minus_f)];
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Power [kW]')
axis([0 T+1 min([min(P_grid_plus_f-P_grid_minus_f),sum(P_grid_plus_f-P_grid_minus_f)])-100 max([max(P_grid_plus_f-P_grid_minus_f),sum(P_grid_plus_f-P_grid_minus_f)])+100])
xticks(0:2:T+1)
Leg = [Leg,"TOTAL"];
legend(Leg,'Location','Best')
title('Main Grid Exchange')

% Wastewater Output
figure(6)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;W_WW_f(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min(min(W_WW_f))-100 max(max(W_WW_f))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeast')
title('Wastewater Output')

% Water Treatment Output
figure(7)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;W_WT_f(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min(min(W_WT_f))-100 max(max(W_WT_f))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Best')
title('Water Treatment Units')

% Water Storage Input/Output
figure(8)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;W_STc_f(m,:)-W_STd_f(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min(min(W_STc_f-W_STd_f))-100 max(max(W_STc_f-W_STd_f))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeast')
title('Water Storage Input/Output')

% Water Storage Level
figure(9)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;WL_ST_f(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water [gal]')
axis([0 T+1 min(min(WL_ST_f))-1000 max(max(WL_ST_f))+1000])
xticks(0:2:T+1)
legend(Leg,'Location','Northeast')
title('Water Storage Levels')

% Main Water System
figure(10)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;W_main_plus_f(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
Y = [Y;sum(W_main_plus_f)];
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min([min(W_main_plus),sum(W_main_plus_f)])-100 max([max(W_main_plus),sum(W_main_plus_f)])+100])
xticks(0:2:T+1)
Leg = [Leg,"TOTAL"];
legend(Leg,'Location','Best')
title('Main Water Import')

% Water Load
figure(11)
Y = [];
Leg = [];
W_ES_total = sum(W_ES_f,3);
for m = 1:M
    Y = [Y;W_LOAD1(m,:)+W_ES_total(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min(min(W_LOAD1+W_ES_total))-100 max(max(W_LOAD1+W_ES_total))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeast')
title('Water Net Demand')
