function [ V, Fx, Cvx ] = fn_static_lateral_force( w, h, T_elastic, Sds, Sd1, model )
% Static Lateral Force Method as outlined in ASCE 41-17

%% Assumptions
% 1) T-long = 8 sec
% 2) Assume Dmax site parameters
% 3) Assume RC SMF, where R = 8, use as mu_strength 

%% Initial setup
% Import Packages
import asce_41.fn_c_factors


%% T_empirical
% hard coded for wall, but not used
Ct = 0.02;
h_n = h(end)/12;
beta = 0.75;
T_emp = Ct * h_n ^ beta;

%% Caclulate seismic response coefficient: Cs
haz_class = model.framing_system_1{1}(1:2);
[ c1, c2, c_m ] = fn_c_factors( model.site_class, model.num_stories, T_elastic, haz_class, [], model.R );

%% Determin Sa from Design Spectra
Sa = min(Sds, Sd1/T_elastic);

%% Calculate Base Shear: V
W = sum(w); % Total effective seismic weight
V = c1*c2*c_m*Sa*W; % Total Base shear eq 7-21

%% Calculate vertical distribution factor: Cvx
k = min([max([interp1([0.5,2.5],[1,2],T_elastic,'linear','extrap'),1]),2]);
Cvx = (w .* h.^k)/sum(w .* h.^k); % eq 7-25

%% Calculate Vertical Forces
Fx = Cvx*V; % eq 7-24

end

