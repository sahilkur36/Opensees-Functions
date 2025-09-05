function [ analysis ] = fn_analysis_options( analysis )
% Description: Defaults secondary analysis options for ASCE 41 Assessment. 

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs: Analysis Data Structure

% Outputs: Analysis Data Structure

% Assumptions:

%% Basic Defaults                    
% Run Options
analysis.run_opensees = 1; % 1 = Run opensees, 0 = use existing results
analysis.run_opensees_post_process = 1; % 1 = Run opensees, 0 = use existing results
analysis.asce_41_post_process = 0; % 1 = run asce 41 post process logic
analysis.opensees_SP = 0; % 0 = Standard OpenSees; 1 = OpenseesSP
analysis.summit = 0; % Write tcl files to be run on summit and change location of opensees call
analysis.skip_2_outputs = 0; % Skip all the way to the plotters

% Model Options
analysis.stories_nonlinear = inf; % Default to all modeling all stories as nonlinear when doing NDP
analysis.stories_nonlinear_low = 0; % all stories at or below this story to be elastic (0 = all nonlinear)
analysis.elastic_beams = 0; % 0 = beams can be nonlinear (default), 1 = beams are assumed to be elastic
analysis.rigid_diaphram = 1; % Default the model to assume rigid diaphrams (0 = non-rigid assuption)
analysis.fiber_walls = 0; % 0 = Lumped plasticity walls, 1 = Fiber walls, 2 = MVLEM
analysis.additional_elements = 1; % this is the leaning column

% Opensees Analysis Options
analysis.ground_motion_scale_factor = 1.0; % Scale the GM amplitude   
analysis.damping = 'rayleigh'; % rayleigh, modal, or simple
analysis.hinge_stiff_mod = 10; % Scale up stiffnes of hinges for a lumped plasticiy model. n value from Ibarra paper.
analysis.run_eigen = 1; % Run the eignen anlayis to get mode shapes and periods for the opensees analysis
analysis.initial_timestep_factor = 1; % reduction from eq timestep to analysis timestep
analysis.solution_algorithm = 1; % Run the opensees solution algorthm which will try different things 
analysis.collapse_drift = 0.10; % stop analysis at this drift and say collapse
analysis.joint_model = 1; % 0 = centerline (almost), 1 = ASCE 41 implicit model, 2 = joint 3D, 3 = rigid (via beam/columns)
analysis.joint_explicit = 0; % 0 = rigid, 1 = model joint nonlinearity (could automate this based on first assessment of joints)
analysis.write_xml = 1; % Write and read opensees out files as xml files (0 = .txt files, which is currently broken, on purpose)
analysis.pushover_num_steps = 200; % Number of steps a pushover will take to get to the dirft limit
analysis.cyclic_pushover_peak_drifts = [0.4, 0.5, 0.6]; % Percent of the final Pushover drift of each cycle
analysis.hinge_group_length = 10;
analysis.filter_accel = 0;
analysis.filter_freq_range = [0.5, 1.5];
if ~isfield(analysis,'simple_recorders')
    analysis.simple_recorders = 0;
end
analysis.algorithm = 'KrylovNewton';
analysis.integrator = 'Newmark 0.5 0.25';

% Visuals and Graphics
analysis.element_plots = 0; % Plot hinge backnones and other per element visualizations
analysis.plot_recordings = 0; % Plot analysis results v recorded results
analysis.play_movie = 1; % Have opensees display a real time graphic of the building and analysis
analysis.movie_scale = 10; % Visual scale of the movie playback
analysis.hinge_stories_2_plot = 3;
analysis.suppress_outputs = 0;

%% Define Proceedure Options
if strcmp(analysis.proceedure,'test')
    analysis.type_list = [1];
    analysis.nonlinear_list = [1];
    analysis.dead_load_list = [1];
    analysis.live_load_list = [1];
    analysis.case_list = {'backbones'};
    analysis.pushover_drift_list_x = [0.01];
    analysis.pushover_drift_list_z = [0.05];
    analysis.accidental_torsion_list = [0];
    analysis.damp_ratio_list = [0.03]; % Analysis damping ratio
    
elseif strcmp(analysis.proceedure,'torsion')
    analysis.type_list = [1, 1];
    analysis.nonlinear_list = [0, 0];
    analysis.dead_load_list = [1, 1];
    analysis.live_load_list = [1, 1];
    analysis.case_list = {'NA', 'torsion_check'};
    analysis.pushover_drift_list_x = [NaN, NaN];
    analysis.pushover_drift_list_z = [NaN, NaN];
    analysis.accidental_torsion_list = [0, 1];
    analysis.damp_ratio_list = [0.028, 0.028];
    
elseif strcmp(analysis.proceedure,'Pushover')
    analysis.type_list = [2]; % NonLinear Pushover then NL Pushover x 2
    analysis.nonlinear_list = [1];
    analysis.dead_load_list = [1];
    analysis.live_load_list = [0.2];
    analysis.case_list = {'backbones'};
    analysis.pushover_drift_list_x = [0.10]; % Drift limit where the pushover will go till
    analysis.pushover_drift_list_z = [0.10];
    analysis.accidental_torsion_list = [0];
    analysis.damp_ratio_list = [0.03]; % Analysis damping ratio
    
elseif strcmp(analysis.proceedure,'Pushover2')
    analysis.type_list = [2, 2]; % NonLinear Pushover then NL Pushover x 2
    analysis.nonlinear_list = [1, 1];
    analysis.dead_load_list = [1, 1];
    analysis.live_load_list = [1, 1];
    analysis.case_list = {'backbones_pushover', 'backbones'};
    analysis.pushover_drift_list_x = [0.05, 0.05]; % Drift limit where the pushover will go till
    analysis.pushover_drift_list_z = [0.05 0.05];
    analysis.accidental_torsion_list = [0, 0];
    analysis.damp_ratio_list = [0.03, 0.03]; % Analysis damping ratio
    
elseif strcmp(analysis.proceedure,'NDP')
    analysis.type_list = [2, 2, 1]; % Linear Pushover then NL Pushover x 2 then 1 NL dynamic
    analysis.nonlinear_list = [1, 1, 1];
    analysis.dead_load_list = [1, 1, 1];
    analysis.live_load_list = [1, 1, 1];
    analysis.case_list = {'backbones_pushover', 'backbones', 'NA'};
    analysis.pushover_drift_list_x = [0.03, 0.03, NaN]; % Drift limit where the pushover will go till
    analysis.pushover_drift_list_z = [0.01, 0.01, NaN];
    analysis.accidental_torsion_list = [0, 0, 1];
    analysis.damp_ratio_list = [0.03, 0.03, 0.03]; % Analysis damping ratio
   
elseif strcmp(analysis.proceedure,'LDP') % Linear Test
    analysis.type_list = [2, 2, 1, 1, 1]; % 1 = dynamic, 2 = pushover % 3 = static cyclic
    analysis.nonlinear_list = [1, 1, 0, 0, 0]; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
    analysis.dead_load_list = [1, 1, 1.1, 0.9, 1.1]; % Dead load factor for linear analysis
    analysis.live_load_list = [1, 1, 1.1, 0.0, 1.1]; % Live load factor for linear analysis
    analysis.case_list = {'backbones_pushover', 'backbones', 'load_case_1', 'load_case_2', 'load_case_3'};
    analysis.pushover_drift_list_x = [0.02, 0.02, NaN, NaN, NaN]; % Drift limit where the pushover will go till
    analysis.pushover_drift_list_z = [0.01, 0.01, NaN, NaN, NaN];
    analysis.accidental_torsion_list = [0, 0, 1, 1, 0];
    analysis.damp_ratio_list = [0.05, 0.05, 0.05, 0.05, 0.05]; % Analysis damping ratio
    
elseif strcmp(analysis.proceedure,'LSP') % Linear Test
    analysis.type_list =          [4, 4]; % 1 = dynamic, 2 = pushover % 3 = static cyclic, 4 =  static loads
    analysis.nonlinear_list =     [0, 0]; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
    analysis.drift_run_list =     [0, 0]; % run exceptions for ELFP drift cases
    analysis.dead_load_list =     [1.1,  0.9]; % Dead load factor for linear analysis
    analysis.live_load_list =     [1.1,  0]; % Live load factor for linear analysis
%     analysis.live_out_load_list = [1.1,  0]; % Live load factor on outer bays for linear analysis
%     analysis.live_in_load_list =  [1.1,  0]; % Live load factor on outer bays for linear analysis
%     analysis.eq_lat_load_list =   [1, 1]; % earthquake load factor for linear analysis --- The 5 is the R factor for LSL assessment for walls
%     analysis.eq_vert_load_list =  [0, 0]; % earthquake vertical load factor for linear analysis
    analysis.case_list = {'load_case_1', 'load_case_2'};
    analysis.accidental_torsion_list = [0, 0];
    analysis.damp_ratio_list = [0.05, 0.05]; % Analysis damping ratio
    
elseif (strcmp(analysis.proceedure,'ELFP') || strcmp(analysis.proceedure,'MRSA')) && ~analysis.run_drifts % Equvalent Lateral Force Procedure - For forces
    analysis.type_list =          [4, 4, 4, 4]; % 1 = dynamic, 2 = pushover % 3 = static cyclic, 4 =  static loads
    analysis.nonlinear_list =     [0, 0, 0, 0]; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
    analysis.drift_run_list =     [0, 0, 0, 0]; % run exceptions for ELFP drift cases
    analysis.dead_load_list =     [1.4, 1.2, 1.2,  0.9]; % Dead load factor for linear analysis
    analysis.live_out_load_list = [0,   1.6, 0.5,  0]; % Live load factor on outer bays for linear analysis
    analysis.live_in_load_list =  [0,   1.6, 0.5,  0]; % Live load factor on outer bays for linear analysis
    analysis.eq_lat_load_list =   [0,   0,   1,    1]; % earthquake load factor for linear analysis --- The 5 is the R factor for LSL assessment for walls
    analysis.eq_vert_load_list =  [0,   0,   0.2, -0.2]; % earthquake vertical load factor for linear analysis
    analysis.case_list = {'load_case_1', 'load_case_2', 'load_case_3', 'load_case_4'};
    analysis.accidental_torsion_list = [0, 0, 0, 0];
    analysis.damp_ratio_list = [0.05, 0.05, 0.05, 0.05]; % Analysis damping ratio
    
elseif (strcmp(analysis.proceedure,'ELFP') || strcmp(analysis.proceedure,'MRSA')) && analysis.run_drifts % Equvalent Lateral Force Procedure - For drifts
    analysis.type_list =          [4]; % 1 = dynamic, 2 = pushover % 3 = static cyclic, 4 =  static loads
    analysis.nonlinear_list =     [0]; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
    analysis.drift_run_list =     [0]; % run exceptions for ELFP drift cases
    analysis.dead_load_list =     [1]; % Dead load factor for linear analysis
    analysis.live_out_load_list = [0.4]; % Live load factor on outer bays for linear analysis
    analysis.live_in_load_list =  [0.4]; % Live load factor on outer bays for linear analysis
    analysis.eq_lat_load_list =   [1]; % earthquake load factor for linear analysis --- The 5 is the R factor for LSL assessment for walls
    analysis.eq_vert_load_list =  [0]; % earthquake vertical load factor for linear analysis
    analysis.case_list = {'load_case_1'};
    analysis.accidental_torsion_list = [0];
    analysis.damp_ratio_list = [0.05]; % Analysis damping ratio
    
end


end

