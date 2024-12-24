function [GHI_av,DNI_av,DHI_av] = SolarIR(SOLAR)
% Solar irradiance calculator using hourly NREL TMY3 solar data
 % Outputs are 24 hr average vectors

GHI = [];%[W/m^2]
DNI = [];%[W/m^2]
DHI = [];%[W/m^2]
for c = 1:24:length(SOLAR)
    GHI = [GHI SOLAR(c:c+23,3)];
    DNI = [DNI SOLAR(c:c+23,6)];
    DHI = [DHI SOLAR(c:c+23,9)];
end
GHI_av = mean(GHI,2);% average GHI at every time interval [W/m^2]
DNI_av = mean(DNI,2);% average DNI at every time interval [W/m^2]
DHI_av = mean(DHI,2);% average DHI at every time interval [W/m^2]

end

