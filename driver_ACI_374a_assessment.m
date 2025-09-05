%% Clear the Workspace
clear
close
clc
fclose('all');

%% Description: Run mutliple types of Opensees assessments for ACI study

% Created By: Dustin Cook
% Date Created: 12/06/2023

% Inputs:

% Outputs:

% Assumptions:

%% Set models and assessments to run
% models = {'m108vLSL'	'm108vALT'	'm108v1'	'm204v1'	'm208v1'	'm208v2'	'm208v3'	'm212v1'	'm212v2'};
% models = {'m108v1', 'm108vLSL', 'm108vALT', 'm112v1'};
% models = {'m104v1', 'm108v1', 'm108vALT', 'm112v1'};
models = {'m104v10', 'm104v12'};
% models = {'m112v1'};
% models = {'m104v1', 'm112v1'};
% models = {'m104v1'};
% models = {'m212v1'};
% models = {'m204v1' 'm204v2' 'm204v3' 'm208v1' 'm208v2' 'm208v3' 'm212v1' 'm212v2' 'm212v3'};
% models = {'m208v1' 'm208v2'  'm208v3'};
% models = {'m208v1'};
% models = {'m208v1'	'm208v2'	'm208v3'};
% procedures = {'MRSA', 'ELFP', 'MRSA', 'ELFP'};
% procedure_subname = {'drifts', 'drifts', 'forces', 'forces'};
% run_drifts = [1, 1, 0, 0];
% nonlinear = [0, 0, 0, 0];
% simple_recorders = [0, 0, 0, 0];
% R_override = [0, 0, 1, 1];

% procedures = {'MRSA', 'MRSA'};
% % procedure_subname = {'forces', 'drifts'};
% procedure_subname = {'drifts', 'forces'};
% % run_drifts = [0, 1];
% run_drifts = [1, 0];
% nonlinear = [0, 0];
% simple_recorders = [0, 0];
% R_override = [1, 1];
% Cd_override = [1, 1];

procedures = {'MRSA'};
procedure_subname = {'drifts'};
run_drifts = [1];
nonlinear = [0];
simple_recorders = [0];
R_override = [1];
Cd_override = [1];

% procedures = {'P58', 'P58', 'ELFP', 'ELFP', 'MRSA'};
% procedure_subname = {'NL', 'Linear' 'drifts', 'forces', ''};
% run_drifts = [0 0 1 0 0];
% nonlinear = [1 0 0 0 0];
% simple_recorders = [1 0 0 0 0];
analysis_id = '2024_10_29-aci';

%% Load archetype model data
model_data = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
scale_factors = [0.25 0.5 0.75 1.0 1.25 1.5 1.75 2];

%% Run each model through each assessment
if 1 == 1
    for m = 1:length(models)
        model = model_data(strcmp(model_data.id,models{m}),:);
        for p = 1:length(procedures)
            analysis = [];
            analysis.clear_existing_data = 0;
            fprintf('Running Model %s through assessment procedure %s\n', model.name{1}, [procedures{p} '-' procedure_subname{p}])
            
            % User Inputs (Think about changing this to a file read and command line execution)
            analysis.model_id = model.id{1};
            analysis.proceedure = procedures{p}; % ELFP or LSP or LDP or NDP or test
            analysis.run_drifts = run_drifts(p);
            analysis.nonlinear = nonlinear(p);
            analysis.simple_recorders = simple_recorders(p);
            analysis.id = [procedure_subname{p} '-' analysis_id]; % ID of the analysis for it to create its own directory
            

            if strcmp(procedures{p},'P58') % Run NLRHA through gm set for FEMA P-85 assessment
                % Import Packages
                import ida.fn_master_P58_response
                
                % Fixed inputs
                remote_dir = [];
                [ analysis ] = fn_analysis_options_P58( analysis );
                
                % Run Assessment
                tic
                fn_master_P58_response( model, analysis, remote_dir )
                toc
                
            else % ASCE 41 Drivers
                % Import Packages
                import asce_41.main_ASCE_41

                % Fixed inputs
                analysis.R_override = R_override(p);
                analysis.Cd_override = Cd_override(p);
                analysis.model_type = 3; % 1 = SDOF, 2 = MDOF, 3 = Archetype model
%                 analysis.eq_lat_load_factor = 1;
                [ analysis ] = fn_analysis_options( analysis );

                % Run Assessment
                for SF = 1:length(scale_factors)
                    tic
                    main_ASCE_41( analysis,  scale_factors(SF) )
                    toc
                end
            end
            
            fprintf('Assessment Complete!\n')
        end
    end
end

% %% Pull Results
% dem_table = table;
% dem_table.story = (1:12)';
% cap_table = table;
% cap_table.story = (1:12)';
% dcr_table = table;
% dcr_table.story = (1:12)';
% model_summary = table;
% for p = 1:length(procedures)
%     remote_dir = ['C:\Users\dtc2\OneDrive - NIST\NIST\ACI 374A\Archetype Model Results' filesep procedures{p} '-' procedure_subname{p}];
%     
%     for m = 1:length(models)
%         model = model_data(strcmp(model_data.id,models{m}),:);
%         
%         % Read tables
%         out_dir = ['outputs' filesep  model.name{1} filesep procedures{p} '_' procedure_subname{p} '-' analysis_id];
%         story = readtable([out_dir filesep 'story.csv']);
%         model_out = readtable([out_dir filesep 'model_data' filesep 'model.csv']);
%         story_model = readtable([out_dir filesep 'model_data' filesep 'story.csv']);
%         story.seismic_wt = story_model.seismic_wt;
% %         story.lateral_force = story_model.lateral_force;
%         
%         % Tabulate DCRs
%         if run_drifts(p)
%             capacity = story.lsl_drift;
%             demand = story.max_drift;
%         else
%             capacity = story.M_lsl;
%             demand = story.Mud;
%         end
%         dem_table.(model.id{1}) = nan(height(dem_table),1);
%         dem_table.(model.id{1})(1:height(story)) = demand;
%         cap_table.(model.id{1}) = nan(height(cap_table),1);
%         cap_table.(model.id{1})(1:height(story)) = capacity;
%         dcr_table.(model.id{1}) = nan(height(dcr_table),1);
%         dcr_table.(model.id{1})(1:height(story)) = demand./capacity;
%         
%         % Tabulate model summary
%         model_summary.model_id{m,1} = model.id{1};
%         model_summary.model_name{m,1} = model.name{1};
%         model_summary.num_stories(m,1) = model.num_stories;
%         model_summary.period(m,1) = model_out.T1_x;
%         model_summary.CuTa(m,1) = model_out.CuTa;
%         model_summary.V(m,1) = model_out.V;
%         model_summary.Cs(m,1) = model_out.Cs;
%         model_summary.max_idr(m,1) = max(story.max_drift);
%         model_summary.max_comp_rot(m,1) = max(story.max_rot);
%         model_summary.max_dcr(m,1) = max(demand./capacity);
%         
%         % plot directory
%         plt_dir = [remote_dir filesep model.id{1}];
%         if ~exist(plt_dir,'dir')
%             mkdir(plt_dir)
%         end
%         
%         % Demand and Capacity Profiles
%         if run_drifts(p)
%             % Plot demand and LSL - idr
%             hold on
%             plot(story.max_drift,story.id)
% %             plot(story.lsl_drift,story.id,'--k')
%             plot(story.lsl_drift*0.8/0.7,story.id,'--k')
%             xlabel('Max Drift')
%             ylabel('Story')
%             box on
%             set(gcf,'position',[0,0,250,400])
%             saveas(gcf,[plt_dir filesep 'idr_profile.png'])
%             close
%             
%             % Plot demand and LSL - rotation
%             hold on
%             plot(story.max_rot,story.id)
%             plot(story.lsl,story.id,'--k')
%             xlabel('Max Rotation')
%             ylabel('Story')
%             box on
%             set(gcf,'position',[0,0,250,400])
%             saveas(gcf,[plt_dir filesep 'rotation_profile.png'])
%             close
%         else
%              % Plot demand and LSL - momemnt
%             hold on
%             plot(story.Mud,story.id)
%             plot(story.M_lsl,story.id,'--k')
%             xlabel('Max Moment')
%             ylabel('Story')
%             box on
%             set(gcf,'position',[0,0,250,400])
%             saveas(gcf,[plt_dir filesep 'moment_profile.png'])
%             close
%         end
%     end
%     
%     % Save data
%     writetable(dem_table,[remote_dir filesep 'dem_table.csv'])
%     writetable(cap_table,[remote_dir filesep 'cap_table.csv'])
%     writetable(dcr_table,[remote_dir filesep 'dcr_table.csv'])
%     writetable(model_summary,[remote_dir filesep 'model_summary.csv'])
%     
% end
