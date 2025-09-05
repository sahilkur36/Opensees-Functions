function [ V, Cs, Fx, Cvx, CuTa ] = fn_equivalent_lateral_force( lat_sys, w, h, T_elastic, S1, Sds, Sd1, R, Ie, run_drifts )
% Equvalent Lateral Force Method as outlined in ASCE 7-16

%% Assumptions
% 1) T-long = 8 sec
% 2) Use Ta for period for forces
% 3) Only applicable to concrete frames and walls currently

%% Intial Parameters
% Import Packages
import asce_7.fn_CuTa

% Fixed Assumptions
T_L = 8;

%% Period Properties
[ CuTa ] = fn_CuTa( lat_sys, Sds, h );
if run_drifts
    T = T_elastic;
else
    T = min(T_elastic,CuTa);
end

%% Caclulate seismic response coefficient: Cs
Cs_short = Sds/(R/Ie); % eq 12.8-2

if T <= T_L
    Cs_T = Sd1/(T*(R/Ie)); % eq 12.8-3
else
    Cs_T = Sd1*T_L/(T^2*(R/Ie)); % eq 12.8-4
end

Cs = min(Cs_short,Cs_T); % Cs need not exceed Cs_T

% eq 12.8-6 ASCE 7-22
if run_drifts
    cs_min_1 = 0; % not requred for drift assessment 
else
    cs_min_1 = max([0.044*Sds*Ie, 0.01]); 
end  

% eq 12.8-7 ASCE 7-22
if S1 >= 0.6
    cs_min_2 = 0.5*S1/(R/Ie);
else
    cs_min_2 = 0;
end

Cs = max([Cs, cs_min_1, cs_min_2]); % eq 12.8-2

%% Calculate Base Shear: V
W = sum(w); % Total effective seismic weight
V = Cs*W; % Total Base shear eq 12.8-1

%% Calculate vertical distribution factor: Cvx
k = min([max([interp1([0.5,2.5],[1,2],T,'linear','extrap'),1]),2]);
Cvx = (w .* h.^k)/sum(w .* h.^k); % eq 12.8-12

%% Calculate Vertical Forces
Fx = Cvx*V; % eq 12.8-11

end

