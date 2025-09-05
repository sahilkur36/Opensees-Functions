function [ ] = fn_master_P58_response( model, analysis, remote_dir )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import Packages
import ida.*
import asce_7.*
import build_model.main_build_model
import opensees.write_tcl.fn_define_model
import opensees.main_eigen_analysis
import usgs.*

% Load basic model data
analysis.model_id = model.id;

% Define read and write directories
main_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id];
tcl_dir = [main_dir '/' 'opensees_data'];
if ~exist(tcl_dir,'dir')
    mkdir(tcl_dir)
end
model_remote_dir = [remote_dir filesep model.id{1}];
if ~exist(model_remote_dir,'dir')
    mkdir(model_remote_dir)
end
if analysis.model_build_remote
    analysis.model_dir = [model_remote_dir '/' 'model_data'];
else
    analysis.model_dir = [main_dir '/' 'model_data'];
    if ~exist(analysis.model_dir,'dir')
        mkdir(analysis.model_dir)
    end
end
pushover_dir = [remote_dir filesep model.name{1}];
%     ELFP_model_dir = ['outputs' '/' model.name{1} '/' 'ELFP' '_' analysis.id '/' 'model_data'];
% tcl_dir = ['outputs' '/' model.name{1} '/' 'ELFP' '_' analysis.id '/' 'eigen_analysis'];
% asce41_dir = [main_dir '/' 'asce_41_data'];

%% Load ground motion data
gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);

%% Set Damping Ratio
%     H = model.story_ht_first + model.story_ht_others*(model.num_stories - 1);
%     analysis.damp_ratio = min(0.36/sqrt(H),0.05); % From equation 3-1 of ATC 114 v1

%% Build Model
analysis.out_dir = main_dir;

% Create model tables
if analysis.model_build_remote == 0
    disp('Building Model ...')
    main_build_model( model, analysis, [] )
end

% Load Model Data
node = readtable([analysis.model_dir filesep 'node.csv'],'readVariableNames',true);
story = readtable([analysis.model_dir filesep 'story.csv'],'readVariableNames',true);
hinge = readtable([analysis.model_dir filesep 'hinge.csv'],'readVariableNames',true);
element = readtable([analysis.model_dir filesep 'element.csv'],'readVariableNames',true);
ele_props_table = readtable([analysis.model_dir filesep 'element_table.csv'],'readVariableNames',true);
joint = readtable([analysis.model_dir filesep 'joint.csv'],'readVariableNames',true);

% Write model tcl file
[ ~ ] = fn_define_model( tcl_dir, node, element, joint, hinge, analysis, model.dimension, story, [], model, ele_props_table );

%% Run Eigen Analysis
[ eigen ] = main_eigen_analysis( analysis, analysis.model_dir );
model.T1_x = eigen.periods(1);

%% Redefine model with fiber
%     analysis.nonlinear_type = 'fiber'; % lumped or fiber
%     main_build_model( model, analysis, [] )
%     node = readtable([model_dir filesep 'node.csv'],'readVariableNames',true);
%     story = readtable([model_dir filesep 'story.csv'],'readVariableNames',true);
%     hinge = readtable([model_dir filesep 'hinge.csv'],'readVariableNames',true);
%     element = readtable([model_dir filesep 'element.csv'],'readVariableNames',true);
%     joint = readtable([model_dir filesep 'joint.csv'],'readVariableNames',true);
%     [ ~ ] = fn_define_model( tcl_dir, node, element, joint, hinge, analysis, model.dimension, story, [], model );

%% Modify Model Data
% Set period variable
ida_results.period = model.T1_x;

% Factor Loads
[ element ] = fn_factor_loads( analysis, element, [] );

%% Set Sa Stripe levels 
analysis.sa_stripes = [];

% Interpolate Site Hazard for building period
if analysis.run_hazards
    [sa_spectra, sa_periods] = fn_call_USGS_hazard_API('E2014', model.lat, model.lng, model.vs30, 1./analysis.rps2run);
    period_cropped = min(max(model.T1_x,min(sa_periods)),max(sa_periods)); % limit building period to the available periods
    sa_levels = interp1(sa_periods,sa_spectra,period_cropped);
    analysis.sa_stripes = [analysis.sa_stripes, sa_levels];
end

% Calc design info
sa_dbe = min(model.sds,model.sd1/model.T1_x);
sa_mce = sa_dbe*1.5;

% design level
if analysis.run_dbe
    analysis.sa_stripes = [analysis.sa_stripes, analysis.de_levels*sa_dbe];
end

% mce level
if analysis.run_mce
    ida_results.mce = sa_mce;
    analysis.sa_stripes = [analysis.sa_stripes, analysis.mce_level*sa_mce];
end

% MAKE SURE THE CONVERSION IS HANGLED PROPERLY ...
%     if period_cropped < 0.2
%         geo_coversion = 1.1;
%     elseif period_cropped < 1
%         geo_coversion = interp1([0.2,1],[1.1,1.3],period_cropped);
%     elseif period_cropped < 5
%         geo_coversion = interp1([1,5],[1.3,1.5],period_cropped);
%     else
%         geo_coversion = 1.5;
%     end
%     sa_mce_geo = sa_mce/geo_coversion;
%     sa_dbe_geo = (2/3)*sa_mce_geo;

%% Run Opensees Models
if analysis.run_ida || analysis.post_process_ida
    fn_master_IDA(analysis, model, story, element, ele_props_table, node, hinge, joint, gm_set_table, ida_results, tcl_dir, main_dir)
end

%% Create Response and Consequence Fragilities
if analysis.create_fragilities
    [col_med, col_beta, p_col_mce, p_col_dbe] = fn_create_collapse_fragility(analysis, gm_set_table, ida_results, main_dir, model_remote_dir, pushover_dir);
    collapse_data.id =  model.id;
    collapse_data.model_name =  model.name;
    collapse_data.med_sa = col_med;
    collapse_data.beta = col_beta;
    collapse_data.p_col_mce = p_col_mce;
    collapse_data.p_col_dbe = p_col_dbe;
end

%% Post Process EDP data
fn_pull_EDPs( main_dir, gm_set_table, story, analysis )

end % funciton

