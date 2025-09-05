%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Assumptions
% 1) 3D model

%% User Inputs
% Define Model
analysis.model_id = 'm104v11';
analysis.proceedure = 'P58'; 
analysis.id = 'NL-2024_12_09'; % ID of the analysis for it to create its own directory
analysis.nonlinear = 1;
analysis.clear_existing_data = 1;

% Define remote directory
% remote_dir = ['C:\Users\dtc2\OneDrive - NIST\NIST\ACI 374A\Models\Archetypes RCMF\EDPs'];
remote_dir = ['C:\Users\dtc2\OneDrive - NIST\NIST\Cost Benefit Study\Models\RCMF - Cook'];

%% Initial Setup
% Import Packages
import ida.fn_master_P58_response

% Define models to run
model_data = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
model_data = model_data(strcmp(model_data.id,analysis.model_id),:);
num_models = height(model_data);

% Set additional analysis parameters
[ analysis ] = fn_analysis_options_P58( analysis );

%% Run and Post Process Analysis
for m = 1:num_models % run for each model     
    fprintf('Running Model %i of %i: %s\n', m, num_models, model_data.name{m})
    model = model_data(m,:);
    fn_master_P58_response( model, analysis, remote_dir )
end
