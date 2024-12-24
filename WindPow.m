function [P_WP] = WindPow(Prated_WT,WS_av,CutIn,v_rated,CutOut)
% Wind power calculated based on wind speed values
%   Calculation of wind power based on wind turbine parameters and
%   metereological wind speed values
    P_WP = [];%[kW]
    for t = 1:1:24
        p_wt = Prated_WT.*((WS_av(t)-CutIn)./(v_rated-CutIn));
        for s1 = 1:1:size(Prated_WT,1)
            for s2 = 1:1:size(Prated_WT,2)
                if WS_av(t)<CutIn(s1,s2) ||WS_av(t)>CutOut(s1,s2)
                    p_wt(s1,s2) = 0;
                elseif WS_av(t)>v_rated(s1,s2) && WS_av(t)<=CutOut(s1,s2)
                    p_wt(s1,s2) = Prated_WT;
                end
            end
        end
        P_WP = [P_WP sum(p_wt,2)];
    end
    
end

