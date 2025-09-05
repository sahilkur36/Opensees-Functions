function [ lateral_force, num_modes ] = fn_MRSA( lat_sys, w, h, eigen, Sds, Sd1, S1, R, Cd, Ie, run_drifts, V_ELFP )
% Modal Response Specturm Method as outlined in ASCE 7-16

%% Assumptions
% 1) T-long = 8 sec
% 2) Use Ta for period for forces
% 3) Only applicable to concrete frames and walls currently

%% Intial Parameters
% Import Packages
import asce_7.fn_CuTa

% Set vars
mode_shapes = eigen.mode_shapes;
periods = eigen.periods;

% Fixed Assumptions
T_L = 8;

%% Period Properties
if ~run_drifts
   [ CuTa ] = fn_CuTa( lat_sys, Sds, h );
   
    % Cap periods at CuTa for forces
    periods(periods > CuTa) = CuTa;
end


%% Get shaking factors
% For constant velocity range
Sa_m = Sd1 ./ periods; % eq 11.4-4

% For long periods
T_L = 8;
filt = periods > T_L;
Sa_m(filt) = Sd1*T_L ./ (periods(filt).^2); % eq 11.4-5

% For short periods
filt = periods <= Sd1/Sds;
Sa_m(filt) = Sds;

% For very short periods...
T0 = 0.2*Sd1/Sds;
filt = periods <= T0;
Sa_m(filt) = Sds*(0.4 + 0.6*periods(filt)/T0); % eq 11.4-3

% Reduce Sa by R
Sa_mD = Sa_m * Ie / R;


Sd = 386*Sa_m.*(periods/(2*pi)).^2;
% Sd = [18.3868792897751	2.52121646901632	0.320450813249867	0.0722657827783561];

%% MRSA according to Lindenburg
Wtot = sum(w);
Lm = sum(w .* mode_shapes);
Mm = sum(w .* mode_shapes .^2);
gamma_n = Lm ./ Mm;
Wm = Lm.^2 ./ Mm;
Vm = Sa_mD .* Wm;

PMm = Wm/Wtot;
for m = 1:length(PMm)
    mass_part = sum(PMm(1:m));
    if mass_part >= 0.9
        num_modes = m;
        break
    end
end

modal_disp = Sd(1:4) .* gamma_n(1:4) .* mode_shapes(:,1:4);

V_srss = sqrt(sum(Vm(:,1:num_modes) .^2,2));
V = max(V_srss,V_ELFP);

if run_drifts
    V = Wtot*0.5*S1/(R/Ie);
    if V >= V_ELFP % Only need to scale drifts if EFLP is controlled by eq 12.8-7
        V_scale_factor = max(V/V_srss,1); % Only scale to eq 12.8-7 of ASCE 7-22
    else
        V_scale_factor = 1; % Otherwise just use MRSA straight
    end
else
    V_scale_factor = max(V_ELFP/V_srss,1); % Scale forces to 100% of ELFP
end

% PMm = Wm / Wtot;
Fx_m = Vm(:,1:num_modes) .* w .* mode_shapes(:,1:num_modes) ./ sum(w .* mode_shapes(:,1:num_modes));
% Fx_srss = sqrt(sum(Fx_m(:,1:mass_idx) .^2,2));
% V_temp = sum(Fx_srss);
% SF = V_ELFP / V_temp;


if run_drifts
    lateral_force = V_scale_factor*Fx_m*(Cd/Ie); % THIS IS WRONG!!!! STRUCTURE SHOULD BE ASSESSED PER MODE... well this is only partially wrong... okay to apply to each mode indivudally before assessing since the scale factor is all the same...
else
    lateral_force = V_scale_factor*Fx_m; % THIS IS WRONG!!!! STRUCTURE SHOULD BE ASSESSED PER MODE... well this is only partially wrong... okay to apply to each mode indivudally before assessing since the scale factor is all the same...
end

% % Fake way... (assuming uniform dist of mass)
% gamma_uniform = sum(mode_shapes) ./ sum(mode_shapes .^2);
% Cx_m = Sa_m .* gamma_uniform .* mode_shapes / model.num_stories;
% Cx_srss = (sum(Cx_m(:,1:mass_idx).^2,2)).^0.5;
% Fx_srss = Cx_srss * Wtot;
% C_srss = sum(Cx_srss);
% V_srss = C_srss * Wtot;

end

