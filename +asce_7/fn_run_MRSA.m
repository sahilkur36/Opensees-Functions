function [ ] = fn_run_MRSA( analysis )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import Packages
import opensees.main_opensees_analysis

% Load model data
model = readtable([analysis.model_dir filesep 'model.csv']);

%% Assess response of each mode
for m = 1:model.num_modes_MRSA
    % Run Opensees Analysis
    if analysis.run_opensees || analysis.run_opensees_post_process
        disp(['Running Mode ' num2str(m)])
        analysis.mode2run = m;
        main_opensees_analysis( analysis )
    end
    
    % Load Analysis Data
    read_dir = [analysis.out_dir filesep 'MRSA' filesep analysis.case filesep 'mode_' num2str(analysis.mode2run) filesep 'opensees_data'];
    load([read_dir filesep 'story_analysis.mat']);
    load([read_dir filesep 'element_analysis.mat']);

    % Save responses from each mode
    story_shear(:,m) = abs(cumsum(story.(['lateral_force_' num2str(m)]),'reverse'));
    max_disp_x(:,m) = story.max_disp_x;
    ave_disp_x(:,m) = story.ave_disp_x;
    Pmax(:,m) = max(element.Pmax,0); % Only compresion forces
    Pmin(:,m) = min(element.Pmin,0); % Only tension forces
    V(:,m) = element.V;
    Mpos(:,m) = element.Mpos;
    Mneg(:,m) = element.Mneg;
    max_rot(:,m) = element.max_rot;
end

% Combine modal responses
story.story_shear = sqrt(sum(story_shear.^2,2));
story.max_disp_x = sqrt(sum(max_disp_x.^2,2));
story.ave_disp_x = sqrt(sum(ave_disp_x.^2,2));
element.Pmax = sqrt(sum(Pmax.^2,2));
element.Pmin = sqrt(sum(Pmin.^2,2));
element.V = sqrt(sum(V.^2,2));
element.M = sqrt(sum(max(abs(Mpos),abs(Mneg)).^2,2));
element.Mpos = sqrt(sum(Mpos.^2,2));
element.Mneg = sqrt(sum(Mneg.^2,2));
element.max_rot = sqrt(sum(max_rot.^2,2));

% Save load case info
save([analysis.out_dir filesep 'MRSA' filesep analysis.case filesep 'story_analysis.mat'],'story')
save([analysis.out_dir filesep 'MRSA' filesep analysis.case filesep 'element_analysis.mat'],'element')
    
% if ~strcmp(analysis.case,'NA')
%     write_dir = [analysis.out_dir filesep analysis.case];
%     fn_make_directory( write_dir )
%     save([write_dir filesep 'story_analysis.mat'],'story')
%     save([write_dir filesep 'element_analysis.mat'],'element')
% end

end

