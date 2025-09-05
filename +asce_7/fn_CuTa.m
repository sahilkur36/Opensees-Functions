function [ CuTa ] = fn_CuTa( lat_sys, Sds, h )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Systme Specific Factors
if strcmp(lat_sys,'C1a')
    % RC Frame Archetypes
    Ct = 0.016; % table 12.8-2
    x = 0.9; % table 12.8-2
    hn = max(h)/12;
elseif strcmp(lat_sys,'C2a')
    % RC Wall Archetypes
    Ct = 0.02; % table 12.8-2
    x = 0.75; % table 12.8-2
    hn = max(h)/12;
else
    error('Lateral System Not Regognized. Update to database driven approach')
end

% Coefficient for upper limit on calculated period table 12.8-1
Cu = interp1([0.1, 0.15, 0.2, 0.3],...
             [1.7,1.6,1.5,1.4],...
             min(max(Sds,0.1),0.3));
         
% Calculate CuTa
Ta = Ct*hn^x;
CuTa = Cu*Ta;
    
end

