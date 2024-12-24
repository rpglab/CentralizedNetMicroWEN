function [P_SP] = SolarPow(DNI,DHI,GHI,n,Lat,S_SP,eff_SP,tilt_a,collect_a)
% Solar power calculator based on irradiance values
%   Calculation of solar power based on solar panel parameters as well as
%   GHI, DNI, and DHI values

    delta = 23.45*sind((360/365)*(n-81));%[°]
    Hour_a = 15*(12-[1:24])';%[°]
    Alt_a = rad2deg(asin(cosd(Lat)*cosd(delta)*cosd(Hour_a) + sind(Lat)*sind(delta)));%[°]
    Solar_a = rad2deg(asin((cosd(delta)*sind(Hour_a))./cosd(Alt_a)));%[°]
    incidence = [];%[rad]
    % Incidence angle- Rows: time intervals ; Columns: solar panels
    for t = 1:1:24
        incidence = [incidence;acos(cosd(Alt_a(t))*cosd(Solar_a(t) - collect_a')*sind(tilt_a) + sind(Alt_a(t))*cosd(tilt_a))];
    end
    incidence = rad2deg(incidence);%[°]
    reflectance = 0.2;% Ground light reflectance [-]
    % Calculating irradiances
    I_BCav = [];%[W/m^2]
    I_DCav = [];%[W/m^2]
    I_RCav = [];%[W/m^2]
    for t = 1:1:24
        I_BCav = [I_BCav;DNI(t)*cosd(incidence(t))];
        I_DCav = [I_DCav;DHI(t)*((1+cosd(tilt_a))/2)];
        I_RCav = [I_RCav;GHI(t)*reflectance*((1-cosd(tilt_a))/2)];
    end
    I_av = (I_BCav + I_DCav + I_RCav)';%[W/m^2]
    % SOLAR POWER OUTPUT
    % rows: buses; columns: time interval
    P_SP = [];%[W]
    for t = 1:1:24
        P_SP = [P_SP, eff_SP*S_SP*I_av(t)];
    end
    P_SP(P_SP<0) = 0;
    P_SP = P_SP/1000;%[kW]

end

