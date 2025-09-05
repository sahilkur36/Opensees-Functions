function [ ] = fn_combine_load_cases( analysis, load_case_id )
% Description: Main script that post process an ELFP runs and combine loads
% from various load combos

% Created By: Dustin Cook
% Date Created: 7/9/2021

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Define Read and Write Directories
if strcmp(analysis.proceedure,'MRSA')
    read_dir = [analysis.out_dir filesep 'MRSA' filesep analysis.case];
else
    read_dir = [analysis.out_dir filesep 'opensees_data'];
end
write_dir = [analysis.out_dir filesep 'asce_7_data'];
if ~exist(write_dir,'dir')
    fn_make_directory( write_dir )
end

% Load Analysis Data
% load([read_dir filesep 'model_analysis.mat'])
load([read_dir filesep 'story_analysis.mat'])
load([read_dir filesep 'element_analysis.mat'])

% % Hard coded var that should be passed in as inputs
% Cd = 5.5;
% Ie = 1;
% 
% %% Caluclate Story Drifts
% disp = Cd * [0; story.ave_disp_x] / Ie;
% story.drift = abs(disp(2:end) - disp(1:(end-1))) ./ story.story_ht;

%% Combine demands from various load cases
% Looks like I am not saving each specific load case here, just overwriting
% each run, and only saving the combined envelope
if load_case_id == 1 % this is the first load case, just save the data
    save([write_dir filesep 'story_analysis.mat'],'story')
    save([write_dir filesep 'element_analysis.mat'],'element')
else % combine with previously run load cases
    % load previously saved demands
    story_combo = load([write_dir filesep 'story_analysis.mat']);
    element_combo = load([write_dir filesep 'element_analysis.mat']);

    % Merge element demands from most recent load case
    element.Pmax = max(abs([element_combo.element.Pmax,element.Pmax]),[],2);
    element.Pmin = max(abs([element_combo.element.Pmin,element.Pmin]),[],2);
    element.V = max(element_combo.element.V,element.V);
%     story_filt = element.story == 2;
%     bm_filt = strcmp(element.type,'beam');
%     element_combo.element.Mpos(story_filt & bm_filt)
%     element.Mpos(story_filt & bm_filt)
    element.Mpos = max(abs([element_combo.element.Mpos,element.Mpos]),[],2);
    element.Mneg = max(abs([element_combo.element.Mneg,element.Mneg]),[],2);

    % Merge story drift demands
    story.ave_disp_x = max(story_combo.story.ave_disp_x,story.ave_disp_x);
    story.max_disp_x = max(story_combo.story.max_disp_x,story.max_disp_x);
    
    % Overwrite saved load case
    save([write_dir filesep 'story_analysis.mat'],'story')
    save([write_dir filesep 'element_analysis.mat'],'element')
end

end
