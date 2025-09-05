function [ ] = main_combine_load_case( analysis, ele_prop_table )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup 
% Import Packages
import asce_41.*

% Load general model data
model = readtable([analysis.model_dir filesep 'model.csv']);

%% Load Various Load Cases
for i = 1:length(analysis.case_list)
    if strcmp(analysis.proceedure,'MRSA')
        read_dir = [analysis.out_dir filesep 'MRSA' filesep analysis.case_list{i}];
    else
        read_dir = [analysis.out_dir filesep analysis.case_list{i}];
    end
    load([read_dir filesep 'story_analysis.mat'])
    load([read_dir filesep 'element_analysis.mat'])
    story_cases.(['load_case_' num2str(i)]) = story;
    element_cases.(['load_case_' num2str(i)]) = element;
end

%% Combine load cases
if strcmp(analysis.proceedure,'LSP') || strcmp(analysis.proceedure,'ELFP') || strcmp(analysis.proceedure,'MRSA') % Combine data for Linear Static Procedure
    
    comp_types = {'column', 'beam', 'wall'};
    comp_props = {'id','story', 'h', 'e', 'a', 'w', 'fc_n', 'iz', 'length'};
%     comp_vars = {'Pmax', 'Pmin', 'V', 'Mpos', 'Mneg', 'max_rot'};
%     comb_method = {'max', 'min', 'max' ,'max', 'min', 'max'};
    comp_vars = {'Pmax', 'Pmin', 'V', 'M', 'max_rot'};
    comb_method = {'max', 'min', 'max' ,'max', 'max'};
    
    % Cobmine element properties
    for c = 1:length(comp_types)
        comp_filt = strcmp(element_cases.load_case_1.type,comp_types{c}) & ~element_cases.load_case_1.rigid;
        prop_filt = ismember(element_cases.load_case_1.Properties.VariableNames,comp_props);
        comp_table.(comp_types{c}) = element_cases.load_case_1(comp_filt,prop_filt);
        for v = 1:length(comp_vars)
            for i = 1:length(analysis.case_list)
                if i == 1
                    comp_table.(comp_types{c}).(comp_vars{v}) = element_cases.(['load_case_' num2str(i)]).(comp_vars{v})(comp_filt);
                else
                    if strcmp(comb_method{v},'max')
                        comp_table.(comp_types{c}).(comp_vars{v}) = max(comp_table.(comp_types{c}).(comp_vars{v}), element_cases.(['load_case_' num2str(i)]).(comp_vars{v})(comp_filt));
                    elseif strcmp(comb_method{v},'min')
                        comp_table.(comp_types{c}).(comp_vars{v}) = min(comp_table.(comp_types{c}).(comp_vars{v}), element_cases.(['load_case_' num2str(i)]).(comp_vars{v})(comp_filt));
                    else
                        error('Unrecognized combination method')
                    end
                end
            end
        end
        
%         % Amplify componnt rotations by Cd
%         if (strcmp(analysis.proceedure,'ELFP') || strcmp(analysis.proceedure,'MRSA')) && analysis.run_drifts
%             comp_table.(comp_types{c}).max_rot = comp_table.(comp_types{c}).max_rot  * model.Cd;
%         end
    end
    
    % Combine Story data
    story_table = table;
    story_table.id = story_cases.load_case_1.id;
    story_table.story_ht = story_cases.load_case_1.story_ht;
    for i = 1:length(analysis.case_list)
        max_drift = (story_cases.(['load_case_' num2str(i)]).max_disp_x  - [0; story_cases.(['load_case_' num2str(i)]).max_disp_x(1:end-1)]) ./ story_table.story_ht;
        if i == 1
            story_table.max_drift = max_drift;
        else
            story_table.max_drift = max(story_table.max_drift, max_drift);
        end
    end
    
%     % Amplify drifts by Cd
%     if (strcmp(analysis.proceedure,'ELFP') || strcmp(analysis.proceedure,'MRSA')) && analysis.run_drifts
%         story_table.max_drift = story_table.max_drift * model.Cd;
%     end
    
    if ~isempty(comp_table.beam) % TEMP ACI 374a LSL calc for beams
        for s = 1:height(story_table)
            bm = comp_table.beam(comp_table.beam.story == story_table.id(s),:);
            story_table.lsl(s) = 0.033 * 0.42;
            story_table.lsl_drift(s) = 0.7 * story_table.lsl(s);
            story_table.max_rot(s) = mean(bm.max_rot);
            story_table.Mud(s) = mean(bm.M);
            story_table.M_lsl(s) = mean(6 .* story_table.lsl(s) .* bm.e .* bm.iz ./ bm.length);
        end
    elseif ~isempty(comp_table.wall) % TEMP ACI 374a LSL calc for walls
        wall = comp_table.wall;
        k1 = 0.1; k2 = 1.2; % for rectangulare walls (Abdullah 2019)
        wall.c_approx = (k1 + k2*(wall.Pmax ./ (wall.a .* wall.fc_n))) .* wall.h;
        story_table.lsl_med = min(max(0.031 - 0.0003 * (wall.h .* wall.c_approx) ./ (wall.w .^2),0.008),0.028);
        story_table.lsl = story_table.lsl_med * 0.42;
        story_table.lsl_drift = 0.7 * story_table.lsl;
        story_table.max_rot = wall.max_rot;
        story_table.Mud = wall.M;
        ls_vec = [];
        for s = 1:height(story)
            ls_vec(s) = (2/3)*sum(story.story_ht(s:end));
        end
        story_table.M_lsl = 3 .* story_table.lsl .* wall.e .* wall.iz ./ ls_vec';
        comp_table.wall = wall;
    else
        error('rework assumptions about systems')
    end
    
    % Save as csv table
    writetable(comp_table.column,[analysis.out_dir filesep 'column_demands.csv'])
    writetable(comp_table.beam,[analysis.out_dir filesep 'beam_demands.csv'])
    writetable(comp_table.wall,[analysis.out_dir filesep 'wall_demands.csv'])
    writetable(story_table,[analysis.out_dir filesep 'story.csv'])
    
    
elseif strcmp(analysis.proceedure,'LDP') % Combine data for Linear Dynamic Procedure
    write_dir = [analysis.out_dir filesep 'asce_41_data'];
    load([write_dir filesep 'node_analysis.mat'])
    load([write_dir filesep 'model_analysis.mat'])

    %% Calculate Envelopes
    % EDP Profiles
    story.max_accel_x = max([story_cases.load_case_1.max_accel_x,story_cases.load_case_2.max_accel_x],[],2);
    story.max_disp_x = max([story_cases.load_case_1.max_disp_x,story_cases.load_case_2.max_disp_x],[],2);
    story.max_disp_center_x = max([story_cases.load_case_1.max_disp_center_x,story_cases.load_case_2.max_disp_center_x],[],2);
    story.max_drift_x = max([story_cases.load_case_1.max_drift_x,story_cases.load_case_2.max_drift_x],[],2);
    if isfield(story,'max_drift_z')
        story.max_accel_z = max([story_cases.load_case_1.max_accel_z,story_cases.load_case_2.max_accel_z],[],2);
        story.max_disp_z = max([story_cases.load_case_1.max_disp_z,story_cases.load_case_2.max_disp_z],[],2);
        story.max_disp_center_z = max([story_cases.load_case_1.max_disp_center_z,story_cases.load_case_2.max_disp_center_z],[],2);
        story.max_drift_z = max([story_cases.load_case_1.max_drift_z,story_cases.load_case_2.max_drift_z],[],2);
    end

    %% Torsional Amplification
    story.A_x = ones(height(story),1);
    story.A_z = ones(height(story),1);

    % Calculate the actual torsional amplification multiplier
    TAR_x = story_cases.load_case_1.torsional_factor_x;
    if sum(strcmp('torsional_factor_z',story_cases.load_case_1.Properties.VariableNames)) > 0
        TAR_z = story_cases.load_case_1.torsional_factor_z;
    end

    % X direction
    if sum(TAR_x > 1.2) > 0
        % For linear analysis calcuate the accidental torional amplification factor
        % (this only applies to forces and displacements caused by accidental
        % torsion)
        for i = 1:length(TAR_x)
            story.A_x(i) = min([max([TAR_x(i)/1.2;1])^2,3]);
        end
    end
    % Amplify Displacements
    story.max_disp_x_mod = max(story.max_disp_x, story.max_disp_x + story_cases.load_case_1.max_disp_x*(max(story.A_x)) - story_cases.load_case_3.max_disp_x*max(story.A_x)); 
    story.max_disp_center_x_mod = max(story.max_disp_center_x, story.max_disp_center_x + story_cases.load_case_1.max_disp_center_x*(max(story.A_x)) - story_cases.load_case_3.max_disp_center_x*max(story.A_x)); 
    story.max_drift_x_mod = max(story.max_drift_x, story.max_drift_x + story_cases.load_case_1.max_drift_x*(max(story.A_x)) - story_cases.load_case_3.max_drift_x*max(story.A_x)); 

    % Z direction
    if exist('TAR_z','var')
        if sum(TAR_z > 1.2) > 0
            % For linear analysis calcuate the accidental torional amplification factor
            % (this only applies to forces and displacements caused by accidental
            % torsion)
            for i = 1:length(TAR_z)
                story.A_z(i) = min([max([TAR_z(i)/1.2;1])^2,3]);
            end
        end
        % Amplify Displacements
        story.max_disp_z_mod = max(story.max_disp_z, story.max_disp_z + story_cases.load_case_1.max_disp_z*(max(story.A_z)) - story_cases.load_case_3.max_disp_z*max(story.A_z)); 
        story.max_disp_center_z_mod = max(story.max_disp_center_z, story.max_disp_center_z + story_cases.load_case_1.max_disp_center_z*(max(story.A_z)) - story_cases.load_case_3.max_disp_center_z*max(story.A_z)); 
        story.max_drift_z_mod = max(story.max_drift_z, story.max_drift_z + story_cases.load_case_1.max_drift_z*(max(story.A_z)) - story_cases.load_case_3.max_drift_z*max(story.A_z)); 
    end

    % SRSS for element demands
    story.A_max = max(story.A_x,story.A_z);

    %% Element Forces
    element.Mmax_1_mod = element.Mmax_1 + element_cases.load_case_1.Mmax_1*(max(story.A_max)) - element_cases.load_case_3.Mmax_1*max(story.A_max); 
    element.Mmax_2_mod = element.Mmax_2 + element_cases.load_case_1.Mmax_2*(max(story.A_max)) - element_cases.load_case_3.Mmax_2*max(story.A_max); 
    element.Vmax_1_mod = element.Vmax_1 + element_cases.load_case_1.Vmax_1*(max(story.A_max)) - element_cases.load_case_3.Vmax_1*max(story.A_max); 
    element.Vmax_2_mod = element.Vmax_2 + element_cases.load_case_1.Vmax_2*(max(story.A_max)) - element_cases.load_case_3.Vmax_2*max(story.A_max); 
    element.Pmax_mod = element.Pmax + element_cases.load_case_1.Pmax*(max(story.A_max)) - element_cases.load_case_3.Pmax*max(story.A_max); 

    %% Acceptance
    % Determine amplificaton factors
    [ model, element, story ] = fn_linear_capacity_and_c_factors( model, story, element );

    % Determine acceptance factors
    [ element ] = main_m_factors( ele_prop_table, element );
    [ hinge ] = fn_linear_hinge_accept(element, ele_prop_table, node);

    %% Save Data
    save([write_dir filesep 'model_analysis.mat'],'model')
    save([write_dir filesep 'story_analysis.mat'],'story')
    save([write_dir filesep 'element_analysis.mat'],'element')
    save([write_dir filesep 'hinge_analysis.mat'],'hinge')
end


end

