function [ c1, c2, c_m ] = fn_c_factors( site_class, num_stories, T, haz_class, DCR_max, R )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Import Packages
import asce_41.fn_cm

% Site Class factor per ASCE 41-17 Section 7.4.1.3.1
if strcmp(site_class,'A') || strcmp(site_class,'B')
    a = 130;
elseif strcmp(site_class,'C')
    a = 90;
else
    a = 60;
end

% Calculate Cm
[ c_m ] = fn_cm( num_stories, T, haz_class );

% Calculate C1 and C2 mod factors
if exist('R','var')
    u_strength = R; % If R is passed in, use it as the ductility capacity
else
    u_strength = max([DCR_max*c_m/1.5,1]);
end
c1 = 1 + (u_strength-1)/(a*T^2); % eqn 7-22
c2 = 1 + (1/800)*((u_strength-1)/T)^2; % eqn 7-23


end

