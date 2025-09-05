function [ ] = main_ASCE_41( analysis, SF )
% Description: Main script facilitating an ASCE 41 teir 3 assessment. 

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Import Packages
import build_model.main_build_model
import opensees.main_eigen_analysis
import opensees.main_opensees_analysis
import asce_41.main_ASCE_41_post_process
import asce_41.main_combine_load_case
import asce_41.fn_define_static_loading
import asce_7.fn_run_MRSA

% Pull in database of available models
if analysis.model_type == 1 % SDOF
    model_table = readtable(['inputs' filesep 'sdof_models.csv'],'ReadVariableNames',true);
    ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
elseif analysis.model_type == 2 % MDOF
    model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
    ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
elseif analysis.model_type == 3 % Archetype
    model_table = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
    ele_prop_table = [];
end

% Select Model for the analysis 
model = model_table(strcmp(model_table.id,analysis.model_id),:);

% Create Analysis Directories
analysis.out_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure '_' analysis.id ];
if ~isempty(SF)
    analysis.out_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure '_' analysis.id filesep 'SF_' strrep(num2str(SF),'.','p') ];
end
% analysis.out_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure '_' analysis.id];

analysis.model_dir = [analysis.out_dir filesep 'model_data'];
 if analysis.run_opensees && ~analysis.skip_2_outputs % Don't clear the file if you don't want to run opensees
    fn_make_directory( analysis.out_dir )
    fn_make_directory( analysis.model_dir )
end

%% Begin Assessment
if ~analysis.skip_2_outputs % Don't skip to plotters
    if analysis.run_opensees
        % Run through all the steps of the procedure
        start_idx = 1;
    else
        % Run through just the end
        start_idx = length(analysis.type_list);
    end
    
    for i = start_idx:length(analysis.type_list)
        analysis.type = analysis.type_list(i);
        analysis.nonlinear = analysis.nonlinear_list(i);
        analysis.dead_load = analysis.dead_load_list(i);
        if strcmp(analysis.proceedure,'ELFP') || strcmp(analysis.proceedure,'MRSA')
            analysis.live_out_load = analysis.live_out_load_list(i);
            analysis.live_in_load = analysis.live_in_load_list(i);
            if ~isempty(SF)
                analysis.eq_lat_load_factor = analysis.eq_lat_load_list(i) * SF;
                analysis.eq_vert_load_factor = analysis.eq_vert_load_list(i) * SF;
            else
                analysis.eq_lat_load_factor = analysis.eq_lat_load_list(i);
                analysis.eq_vert_load_factor = analysis.eq_vert_load_list(i);
            end
        else
            analysis.live_load = analysis.live_load_list(i);
        end
        analysis.case = analysis.case_list{i};
        if isfield(analysis,'pushover_drift_list_x')
            analysis.pushover_drift_x = analysis.pushover_drift_list_x(i);
            analysis.pushover_drift_z = analysis.pushover_drift_list_z(i);
        end
        analysis.accidental_torsion = analysis.accidental_torsion_list(i);
        analysis.damp_ratio = analysis.damp_ratio_list(i);
        disp(['Running ' analysis.proceedure ' step ' num2str(i) ' of ' num2str(length(analysis.type_list)) ' ...'])

        %% Build Model
        disp('Building Model ...')
        main_build_model( model, analysis, ele_prop_table )
        
        %% Run Eigen Analysis
        if analysis.run_eigen
            [ eigen ] = main_eigen_analysis( analysis, analysis.model_dir );
        else
            eigen = [];
        end
        
        %% Determine Lateral loads for Static Linear assessments
        if analysis.type ~= 1 % for non-dynamic linear assessments
            fn_define_static_loading( analysis, eigen );
        end
            
        if strcmp(analysis.proceedure,'MRSA') % linear dynamic analysis per ASCE 7
            %% Run each mode of MRSA assessement and combine respones
            fn_run_MRSA( analysis )
        else
            %% Run and Postprocess Opensees Analysis
            disp('Running Opensees ...')
            if analysis.run_opensees || analysis.run_opensees_post_process
                main_opensees_analysis( analysis )
            end

            %% Postprocess ASCE 41 data
            disp('Post Processing Via ASCE 41 ...')
    %         [ capacity(:,i), torsion{i} ] = main_ASCE_41_post_process( analysis, ele_prop_table );
            main_ASCE_41_post_process( analysis, ele_prop_table );

            %% Analysis Checks
    %         disp('Validating Analysis Results ...')
    %         main_check_analysis( analysis, ele_prop_table, capacity, torsion, i )
        end
    end
end

%% Combine Load Cases
% if strcmp(analysis.proceedure,'LDP') || strcmp(analysis.proceedure,'LSP')
    main_combine_load_case( analysis, ele_prop_table )
% end

%% Compile Results and Create Visuals
disp('Plotting Analysis Results ...')
main_plot_analysis_results( analysis, ele_prop_table )

%% Write Element Tables
if strcmp(analysis.proceedure,'NDP')
    fn_pull_ele_database(analysis, ele_prop_table)
end

%% LaTeX Report Writer


disp('Analysis Complete!')
end

