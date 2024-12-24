%% Networked Micro Water Energy Nexus Optimization Model
% Jesus Silva Rodrigues | University of Houston Department of Electrical
% and Computer Engineering

% NetMicroWEN model for single entity controlled network. There are no
% independent system operators, however, system implements a
% post-processing proportional algorithm to adjust exchange values for fair
% and equitative economic benefits for each system in the network. A
% reliability comparison can be made via the power and water load shed
% variables

clear
close all
clc
%% Datasets
% LOADING DATA
SOLAR_coast = readmatrix('SOLAR_GAL.CSV','Range','C5091:M5834');
SOLAR_city = readmatrix('SOLAR_HOU.CSV','Range','C5091:M5834');
WIND_coast = readmatrix('METEO_GAL.xlsx','Range','G5081:G5824');
WIND_city = readmatrix('METEO_HOU.xlsx','Range','G50333:G57772');
P_LOAD = 1000*readmatrix('P_LOAD.xlsx','Range','B4:Y7');%[kW]
W_LOAD = readmatrix('WATER.xlsx','Range','B4:Y7');%[gal/h]

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

% GENERATORS
% Matrix dimensions: (num of MGs, num of gens on each MG) 
P_G_min = [100 150 0;100 150 260;250 0 0;0 0 0];%[kW]
P_G_max = [650 800 0;650 800 940;900 0 0;0 0 0];%[kW]
C_G = [0.33 0.28 0;0.33 0.28 0.23;0.26 0 0;0 0 0];%[$/kWh]
SU_G = [15.00 13.00 0;15.00 13.00 10.35;11.85 0 0;0 0 0];%[$]
NL_G = [11.00 8.70 0;11.00 8.70 7.40;8.45 0 0;0 0 0];%[$/h]
u_G_init = zeros(size(P_G_min));

% ENERGY STORAGE
% Matrix dimensions: (num of MGs, num of ESS on each MG) 
P_ES_limit = [300 300 300 300 300 300;725 725 725 725 0 0;300 300 300 300 300 0;725 725 725 725 725 0];%[kW]
EL_ES_min = [0 0 0 0 0 0;225 225 225 225 0 0;0 0 0 0 0 0;225 225 225 225 225 0];%[kWh]
EL_ES_max = [1660 1660 1660 1660 1660 1660;2500 2500 2500 2500 0 0;1660 1660 1660 1660 1660 0;2500 2500 2500 2500 2500 0];%[kWh]
EL_ES_init = EL_ES_min;
eff_ES_c = [0.80 0.80 0.80 0.80 0.80 0.80;0.95 0.95 0.95 0.95 1 1;0.80 0.80 0.80 0.80 0.80 1;0.95 0.95 0.95 0.95 0.95 1];%[-]
eff_ES_d = [0.60 0.60 0.60 0.60 0.60 0.60;0.98 0.98 0.98 0.98 1 1;0.60 0.60 0.60 0.60 0.60 1;0.98 0.98 0.98 0.98 0.98 1];%[-]
CW_ES = 0.0707*ones(size(P_ES_limit));%[gal/kWh]
CW_ES(P_ES_limit == 0) = 0;
CW_ES(2,:) = 0;% MWEN 2 and MWEN 4 have BESS
CW_ES(4,:) = 0;

% NETWORK POWER TIE-LINE LIMIT
P_E_limit = 1400*ones(4,1);%[kW]

% WASTEWATER TREATMENT
%Reclaimed Wastewater
w_rec = 0.50;% Percentaje of water reclamed
% Wastewater treatment characteristics
W_WW_min = [20;50;20;50];%[gal/h]
W_WW_max = [720;1000;720;1000];%[gal/h]
WL_WW_max = [12000;18600;12000;18600];%[gal]
WL_WW_init = 0.04*WL_WW_max;%[gal]
NL_WW = [75.00;75.00;75.00;75.00];%[$/h]
CW_WW = [365;382;365;382];%[gal/kWh]
eff_WW_pump = 0.84;%[-]

% CONTROLLABLE WATER TREATMENT
W_WT_min = [120;35;150;0];%[gal/h]
W_WT_max = [1500;525;1730;0];%[gal/h]
NL_WT = [17.20;25.5;15.00;0];%[$/h]
CW_WT = [96;4400;100;0];%[gal/kWh]
eff_WT_pump = 0.84;%[-]

% WATER STORAGE
W_ST_limit = [900;900;900;900];%[gal/h]
WL_ST_limit = [10209;10209;10209;10209];%[gal]
WL_ST_init = zeros(size(W_ST_limit));%#[gal]
eff_ST_pump = 0.84;%[-]

% NETWORK WATER TIE-LINE LIMIT
W_E_limit = 980*ones(4,1);%[gal/h]

% SOLAR POWER
S_SP = [10*50;15*180;10*40;15*220];%[m^2]
eff_SP = [0.30;0.35;0.30;0.35];%[-]
tilt_a = 27*ones(size(S_SP));%[°]
collect_a = zeros(size(S_SP));%[°]
Lat_coast = readmatrix('SOLAR_GAL.csv','Range','E1:E1');%[°]
Lat_city = readmatrix('SOLAR_HOU.csv','Range','E1:E1');%[°]
n = 224;% A day in August
P_SP = zeros(size(P_LOAD));
for m = 1:2:size(P_SP,1)
    P_SP(m,:) = SolarPow(DNI_coast,DHI_coast,GHI_coast,n,Lat_coast,S_SP(m),eff_SP(m),tilt_a(m),collect_a(m));
end
for m = 2:2:size(P_SP,1)
    P_SP(m,:) = SolarPow(DNI_city,DHI_city,GHI_city,n,Lat_city,S_SP(m),eff_SP(m),tilt_a(m),collect_a(m));
end

% WIND POWER
Prated_WT = [250 250 300 300;0 0 0 0;250 250 250 250;0 0 0 0];%[kW]
CutIn = [2.4 2.4 3 3;0 0 0 0;2.4 2.4 2.4 2.4;0 0 0 0;];%[m/s]
v_rated = [5.6 5.6 7.2 7.2;0 0 0 0;5.6 5.6 5.6 5.6;0 0 0 0];%[m/s]
CutOut = [10 10 13.6 13.6;0 0 0 0;10 10 10 10;0 0 0 0];%[m/s]
P_WP = zeros(size(P_LOAD));%[kW]
for m = 1:2:size(Prated_WT,1)
    P_WP(m,:) = WindPow(Prated_WT(m,:),WS_coast_av,CutIn(m,:),v_rated(m,:),CutOut(m,:));
end
for m = 2:2:size(Prated_WT,1)
    P_WP(m,:) = WindPow(Prated_WT(m,:),WS_city_av,CutIn(m,:),v_rated(m,:),CutOut(m,:));
end

% GRID AND NETWORK EXCHANGE
% Power
C_grid_plus = Grid';%[$/kWh]
C_grid_minus = 0.8*C_grid_plus;%[$/kWh]
C_Np = 0.9*C_grid_plus;%[$/kWh]
P_grid_limit = (4/4)*sum(P_E_limit);%[kW]
Cp_shed = 1e6;% Power load shed penalty
% Water
C_Wmain_plus = 0.006;%[$/gal]
C_Nw = 0.9*C_Wmain_plus;%[$/gal]
W_main_limit = (4/4)*sum(W_E_limit);%[gal/h]
Cw_shed = 1e6;% Water load shed penalty

%% Centralized Optimization
[M,T] = size(P_LOAD);
G = size(P_G_max,2);
ES = size(P_ES_limit,2);

cvx_begin
% VARIABLES
    % Power Variables
    variable u_G(M,T,G) binary
    variable v_G(M,T,G) binary
    variables P_G(M,T,G) P_grid_plus(M,T) P_grid_minus(M,T)
    variables P_ESc(M,T,ES) P_ESd(M,T,ES) EL_ES(M,T,ES) W_ES(M,T,ES)
    variable e_ESc(M,T,ES) binary
    variable e_ESd(M,T,ES) binary
    variable p_plus(1,T) binary
    variable p_minus(1,T) binary
    variable pn_plus(M,T) binary
    variable pn_minus(M,T) binary
    variables P_N(M,T) P_E(M,T) P_shed(M,T)
    % Water Variables
    variables W_WW(M,T) W_WT(M,T) P_WT(M,T) P_WW(M,T) WL_WW(M,T) WR_WW(M,T)
    variable u_WT(M,T) binary
    variable u_WW(M,T) binary
    variables W_STc(M,T) W_STd(M,T) WL_ST(M,T)
    variable sp_ST(M,T) binary
    variable sv_ST(M,T) binary
    variable an_plus(M,T) binary
    variable an_minus(M,T) binary
    variables W_main_plus(M,T) P_WT_pump(M,T) P_ST_pump(M,T) P_WW_pump(M,T)
    variables W_N(M,T) W_E(M,T) W_shed(M,T)
    
% OBJECTIVE FUNCTION
    % Power cost
    f_E = 0;
    for m = 1:1:M
        for t = 1:1:T
            for g = 1:1:G
                f_E = f_E + SU_G(m,g)*v_G(m,t,g) + delta_t*(NL_G(m,g)*u_G(m,t,g)+C_G(m,g)*P_G(m,t,g)); 
            end
            f_E = f_E + C_grid_plus(t)*P_grid_plus(m,t) - C_grid_minus(t)*P_grid_minus(m,t) + C_Np(t)*P_N(m,t) + Cp_shed*P_shed(m,t);
        end
    end
    % Water cost
    f_W = 0;
    for m = 1:1:M
        for t = 1:1:T
            f_W = f_W + delta_t*(NL_WW(m)*u_WW(m,t)+NL_WT(m)*u_WT(m,t)) + C_Wmain_plus*W_main_plus(m,t) + C_Nw*W_N(m,t) + Cw_shed*W_shed(m,t);
        end
    end
    minimize f_E + f_W
    
% CONSTRAINTS
    subject to
        % Power Constraints
        for t = 1:1:T
            for m = 1:1:M
                % Generators
                for g = 1:1:G
                    P_G_min(m,g)*u_G(m,t,g) <= P_G(m,t,g) <= P_G_max(m,g)*u_G(m,t,g)
                    if t == 1
                        v_G(m,t,g) >= u_G(m,t,g) - u_G_init(m,g)
                    else
                        v_G(m,t,g) >= u_G(m,t,g) - u_G(m,t-1,g)
                    end
                end
                % Energy Storage
                for b = 1:1:size(P_ES_limit,2)
                    0 <= P_ESc(m,t,b) <= P_ES_limit(m,b)*e_ESc(m,t,b)
                    0 <= P_ESd(m,t,b) <= P_ES_limit(m,b)*e_ESd(m,t,b)
                    e_ESc(m,t,b) + e_ESd(m,t,b) <= 1
                    if t == 1
                        EL_ES(m,t,b) == EL_ES_init(m,b) + delta_t*(eff_ES_c(m,b)*P_ESc(m,t,b) - P_ESd(m,t,b)/eff_ES_d(m,b))
                    else
                        EL_ES(m,t,b) == EL_ES(m,t-1,b) + delta_t*(eff_ES_c(m,b)*P_ESc(m,t,b) - P_ESd(m,t,b)/eff_ES_d(m,b))
                    end
                    EL_ES_min(m,b) <= EL_ES(m,t,b) <= EL_ES_max(m,b)
                    W_ES(m,t,b) == CW_ES(m,b)*P_ESc(m,t,b)
                end
                % MicroWEN Power Balance
                sum(P_G(m,t,:),3) + sum(P_ESd(m,t,:)-P_ESc(m,t,:),3) + P_grid_plus(m,t) - P_grid_minus(m,t) + P_N(m,t) == P_LOAD(m,t) - P_SP(m,t) - P_WP(m,t) + P_WW(m,t) + P_WW_pump(m,t) + P_WT(m,t) + P_WT_pump(m,t) + P_ST_pump(m,t) - P_shed(m,t)
                % Power Exchange
                P_E(m,t) == P_grid_plus(m,t) - P_grid_minus(m,t) + P_N(m,t)
                -P_E_limit(m) <= P_E(m,t) <= P_E_limit(m)
                -P_E_limit(m)*pn_minus(m,t) <= P_N(m,t) <= P_E_limit(m)*pn_plus(m,t)
                pn_plus(m,t) + pn_minus(m,t) <= 1
                % Main Grid
                0 <= P_grid_plus(m,t) <= P_E_limit(m)*pn_plus(m,t)
                0 <= P_grid_minus(m,t) <= P_E_limit(m)*pn_minus(m,t)
            end
            0 <= sum(P_grid_plus(:,t)) <= P_grid_limit*p_plus(t)
            0 <= sum(P_grid_minus(:,t)) <= P_grid_limit*p_minus(t)
            p_plus(t) + p_minus(t) <= 1
            % Network Power Balance
            sum(P_N(:,t)) == 0
        end
        % Water Constraints
        for t = 1:1:T
            for m = 1:1:M
                % Wastewater Treatment
                W_WW_min(m)*u_WW(m,t) <= W_WW(m,t) <= W_WW_max(m)*u_WW(m,t)
                if t == 1
                    WL_WW(m,t) == WL_WW_init(m) + delta_t*(WR_WW(m,t)-W_WW(m,t))
                else
                    WL_WW(m,t) == WL_WW(m,t-1) + delta_t*(WR_WW(m,t)-W_WW(m,t))
                end
                0 <= WL_WW(m,t) <= WL_WW_max(m)
                if t == 1
                    WR_WW(m,t) == w_rec*(W_LOAD(m,end))
                else
                    WR_WW(m,t) == w_rec*(W_LOAD(m,t-1)-W_shed(m,t-1))
                end
                W_WW(m,t) == CW_WW(m)*P_WW(m,t)
                % Regular Water Treatment
                W_WT_min(m)*u_WT(m,t) <= W_WT(m,t) <= W_WT_max(m)*u_WT(m,t)
                W_WT(m,t) == CW_WT(m)*P_WT(m,t)
                P_WT(m,t) >= 0
                % Water Storage
                0 <= W_STc(m,t) <= W_ST_limit(m)*sp_ST(m,t)
                0 <= W_STd(m,t) <= W_ST_limit(m)*sv_ST(m,t)
                sp_ST(m,t) + sv_ST(m,t) <= 1
                if t == 1
                    WL_ST(m,t) == WL_ST_init(m) + delta_t*(W_STc(m,t)-W_STd(m,t))
                else
                    WL_ST(m,t) == WL_ST(m,t-1) + delta_t*(W_STc(m,t)-W_STd(m,t))
                end
                0 <= WL_ST(m,t) <= WL_ST_limit(m)
                % MicroWEN Water Balance
                W_WT(m,t) + W_STd(m,t) - W_STc(m,t) + W_main_plus(m,t) + W_N(m,t) == W_LOAD(m,t) - W_WW(m,t) + sum(W_ES(m,t,:)) - W_shed(m,t)
                % Water Exchange
                W_E(m,t) == W_main_plus(m,t) + W_N(m,t)
                -W_E_limit(m) <= W_E(m,t) <= W_E_limit(m)
                -W_E_limit(m)*an_minus(m,t) <= W_N(m,t) <= W_E_limit(m)*an_plus(m,t)
                an_plus(m,t) + an_minus(m,t) <= 1
                % Wastewater Pump
                1000*eff_WW_pump*P_WW_pump(m,t) == 2.3177*W_WW(m,t)
                % Water Treatment Pump
                1000*eff_WT_pump*P_WT_pump(m,t) == 2.3177*W_WT(m,t)
                % Water Storage Pump
                1000*eff_ST_pump*P_ST_pump(m,t) == 2.3177*W_STc(m,t)
                % Main Water System
                0 <= W_main_plus(m,t) <= W_E_limit(m)*an_plus(m,t) 
            end
            0 <= sum(W_main_plus(:,t)) <= W_main_limit
            % Network Water Balance
            sum(W_N(:,t)) == 0  
        end
        % Energy and Water Storage Reliability
        for m = 1:M
            sum(EL_ES(m,T,:)) == 5*0.3*max(P_LOAD(m,:))
            WL_ST(m,T) == 5*0.3*max(W_LOAD(m,:))
            WL_WW(m,T) <= 0.04*WL_WW_max(m)
        end
        % Power and Water Load Shed
        P_shed >= 0
        W_shed >= 0
        
cvx_end

%% Post-Processing Proportional Exchange Adjustment (PEA)

% Power exchange adjustment
P_E_plus = zeros(M,T);
P_E_minus = zeros(M,T);
for t = 1:1:T
    for m = 1:1:M
        if P_E(m,t) > 0
            P_E_plus(m,t) = P_E(m,t);
            P_E_minus(m,t) = 0;
        else
            P_E_minus(m,t) = P_E(m,t);
            P_E_plus(m,t) = 0;
        end
    end
    for m = 1:1:M
        if P_E(m,t) > 0
            if abs(sum(P_E_plus(:,t))) > abs(sum(P_E_minus(:,t)))
                P_N(m,t) = (P_E_plus(m,t)/abs(sum(P_E_plus(:,t))))*abs(sum(P_E_minus(:,t)));
            else
                P_N(m,t) = P_E(m,t);
            end
            P_grid_plus(m,t) = P_E(m,t) - P_N(m,t);
            P_grid_minus(m,t) = 0;
        else
            if abs(sum(P_E_plus(:,t))) < abs(sum(P_E_minus(:,t)))
                P_N(m,t) = (P_E_minus(m,t)/abs(sum(P_E_minus(:,t))))*abs(sum(P_E_plus(:,t)));
            else
                P_N(m,t) = P_E(m,t);
            end
            P_grid_minus(m,t) = abs(P_E(m,t) - P_N(m,t));
            P_grid_plus(m,t) = 0;
        end
    end         
end

% Water exchange adjustment
W_E_plus = zeros(M,T);
W_E_minus = zeros(M,T);
for t = 1:1:T
    for m = 1:1:M
        if W_E(m,t) > 0
            W_E_plus(m,t) = W_E(m,t);
            W_E_minus(m,t) = 0;
        else
            W_E_minus(m,t) = W_E(m,t);
            W_E_plus(m,t) = 0;
        end
    end
    for m = 1:1:M
        if W_E(m,t) > 0
            if abs(sum(W_E_plus(:,t))) > abs(sum(W_E_minus(:,t)))
                W_N(m,t) = (W_E_plus(m,t)/abs(sum(W_E_plus(:,t))))*abs(sum(W_E_minus(:,t)));
            else
                W_N(m,t) = W_E(m,t);
            end
            W_main_plus(m,t) = W_E(m,t) - W_N(m,t);
        else
            if abs(sum(W_E_plus(:,t))) < abs(sum(W_E_minus(:,t)))
                W_N(m,t) = (W_E_minus(m,t)/abs(sum(W_E_minus(:,t))))*abs(sum(W_E_plus(:,t)));
            else
                W_N(m,t) = W_E(m,t);
            end
            W_main_plus(m,t) = 0;
        end
    end         
end

%% Updated Objective Values

% Power cost
f_E = zeros(M,1);
for m = 1:1:M
    for t = 1:1:T
        for g = 1:1:G
            f_E(m) = f_E(m) + SU_G(m,g)*v_G(m,t,g) + delta_t*(NL_G(m,g)*u_G(m,t,g)+C_G(m,g)*P_G(m,t,g));
        end
        f_E(m) = f_E(m) + C_grid_plus(t)*P_grid_plus(m,t) - C_grid_minus(t)*P_grid_minus(m,t) + C_Np(t)*P_N(m,t);
    end
end

% Water cost
f_W = zeros(M,1);
for m = 1:1:M
    for t = 1:1:T
        f_W(m) = f_W(m) + delta_t*(NL_WW(m)*u_WW(m,t)+NL_WT(m)*u_WT(m,t)) + C_Wmain_plus*W_main_plus(m,t) + C_Nw*W_N(m,t);
    end
end

% Total cost
f_cost = f_E + f_W;

%% Output Results
fprintf("Optimal Operation Costs\n")
for m = 1:M
    fprintf("MWEN%d: %0.2f\n",m,f_cost(m))
end
fprintf("TOTAL: %0.2f.\n",sum(f_cost))

% Net Loads
figure(1)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;P_LOAD(m,:)-P_SP(m,:)-P_WP(m,:)+P_WW(m,:)+P_WW_pump(m,:)+P_WT(m,:)+P_WT_pump(m,:)+P_ST_pump(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Power [kW]')
axis([0 T+1 min(min(P_LOAD-P_WP-P_SP+P_WW+P_WW_pump+P_WT+P_WT_pump+P_ST_pump))-100 max(max(P_LOAD-P_WP-P_SP+P_WW+P_WW_pump+P_WT+P_WT_pump+P_ST_pump))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeastoutside')
title('Power Net Demand')

% Energy Storage charging/discharging
figure(2)
Y = [];
Leg = [];
P_ESc_all = sum(P_ESc,3);
P_ESd_all = sum(P_ESd,3);
for m = 1:M
    Y = [Y;P_ESc_all(m,:)-P_ESd_all(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Power [kW]')
legend(Leg,'Location','Northeastoutside')
axis([0 T+1 min(min(P_ESc_all-P_ESd_all))-100 max(max(P_ESc_all-P_ESd_all))+100])
xticks(0:2:T+1)
title('Energy Storage Input/Output')

% Energy Storage charge level
figure(3)
Y = [];
Leg = [];
EL_ES_all = sum(EL_ES,3);
for m = 1:M
    Y = [Y;EL_ES_all(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Charge Level [kWh]')
legend(Leg,'Location','Northeastoutside')
axis([0 T+1 min(min(EL_ES_all))-1000 max(max(EL_ES_all))+1000])
xticks(0:2:T+1)
title('Energy Storage Charge Level')

% Generators Output
figure(4)
Y = [];
Leg = [];
P_G_all = sum(P_G,3);
for m = 1:M
    Y = [Y;P_G_all(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Power [kW]')
legend(Leg,'Location','Northeastoutside')
axis([0 T+1 min(min(P_G_all))-100 max(max(P_G_all))+100])
xticks(0:2:T+1)
title('Generators')

% Main Grid
figure(5)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;P_grid_plus(m,:)-P_grid_minus(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Power [kW]')
axis([0 T+1 min(min(P_grid_plus-P_grid_minus))-100 max(max(P_grid_plus-P_grid_minus))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeastoutside')
title('Main Grid Exchange')

% Network Power Exchange
figure(6)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;P_N(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Power [kW]')
axis([0 T+1 min(min(P_N))-100 max(max(P_N))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Southwest')
%title('Power Echange Within Network')

% Wastewater Output
figure(7)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;W_WW(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min(min(W_WW))-100 max(max(W_WW))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeastoutside')
title('Wastewater Output')

% Water Treatment Output
figure(8)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;W_WT(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min(min(W_WT))-100 max(max(W_WT))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeastoutside')
title('Water Treatment Units')

% Water Storage Input/Output
figure(9)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;W_STc(m,:)-W_STd(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min(min(W_STc-W_STd))-100 max(max(W_STc-W_STd))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeastoutside')
title('Water Storage Input/Output')

% Water Storage Level
figure(10)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;WL_ST(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water [gal]')
axis([0 T+1 min(min(WL_ST))-1000 max(max(WL_ST))+1000])
xticks(0:2:T+1)
legend(Leg,'Location','Northeastoutside')
title('Water Storage Levels')

% Main Water System
figure(11)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;W_main_plus(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min(min(W_main_plus))-100 max(max(W_main_plus))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeastoutside')
title('Main Water Import')

% Network Water Exchange
figure(12)
Y = [];
Leg = [];
for m = 1:M
    Y = [Y;W_N(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min(min(W_N))-100 max(max(W_N))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeastoutside')
%title('Water Echange Within Network')

% Water Load
figure(13)
Y = [];
Leg = [];
W_ES_total = sum(W_ES,3);
for m = 1:M
    Y = [Y;W_LOAD(m,:)+W_ES_total(m,:)];
    Leg = [Leg,sprintf("MWEN_%d",m)];
end
bar(1:T,Y)
grid on
xlabel('Time Interval [h]')
ylabel('Water Rate [gal/h]')
axis([0 T+1 min(min(W_LOAD+W_ES_total))-100 max(max(W_LOAD+W_ES_total))+100])
xticks(0:2:T+1)
legend(Leg,'Location','Northeastoutside')
title('Water Net Demand')
