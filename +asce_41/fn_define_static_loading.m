function [ ] = fn_define_static_loading( analysis, eigen )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Import Packages
import asce_7.*
import asce_41.fn_static_lateral_force

% Load story data
model = readtable([analysis.model_dir filesep 'model.csv']);
story = readtable([analysis.model_dir filesep 'story.csv']);

% equivalent lateral force vertical distribution
w = story.seismic_wt;
h = story.story_ht + story.y_start;
if isfield(analysis,'R_override') && analysis.R_override > 0
    R = analysis.R_override;
else
    R = model.R;
end
if isfield(analysis,'Cd_override') && analysis.Cd_override > 0
    Cd = analysis.Cd_override; 
else
    Cd = model.Cd;
end

if ~isfield(analysis,'run_drifts')
    analysis.run_drifts = 0; %default not to run if it hasn't been defined
end
if strcmp(analysis.proceedure,'LSP') || strcmp(analysis.proceedure,'Pushover')
    [ ~, story.lateral_force, ~ ] = fn_static_lateral_force( w, h, model.T1_x, model.sds, model.sd1, model );
elseif strcmp(analysis.proceedure,'ELFP') 
    [ model.V, model.Cs, story.lateral_force, ~ , model.CuTa ] = fn_equivalent_lateral_force( model.framing_system_1, w, h, model.T1_x, model.s1, model.sds, model.sd1, R, model.ie, analysis.run_drifts );
elseif strcmp(analysis.proceedure,'MRSA')
    % Get target base shear from ELFP
    [ V_ELFP, model.Cs, ~, ~ , model.CuTa] = fn_equivalent_lateral_force( model.framing_system_1, w, h, model.T1_x, model.s1, model.sds, model.sd1, R, model.ie, analysis.run_drifts );
    model.V = V_ELFP;

    % Perform Modal Analysis
    [ lateral_force, model.num_modes_MRSA ] = fn_MRSA( model.framing_system_1, w, h, eigen, model.sds, model.sd1, model.s1, R, Cd, model.ie, analysis.run_drifts, V_ELFP );
    for n = 1:model.num_modes_MRSA
        story.(['lateral_force_' num2str(n)]) = lateral_force(:,n);
    end
end

%% Save Select Data
writetable(story,[analysis.model_dir filesep 'story.csv'])
writetable(model,[analysis.model_dir filesep 'model.csv'])

end

