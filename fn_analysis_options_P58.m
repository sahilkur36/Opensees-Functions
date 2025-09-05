function [ analysis ] = fn_analysis_options_P58( analysis )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Assessment Parameters
analysis.model_type = 3; % 1 = SDOF, 2 = MDOF (default), 3 = Archetype model
analysis.type = 1; % dynamic analysis
analysis.dead_load = 1; % Dead load factor
analysis.live_load = 0.25; % live load factor
analysis.collapse_drift = 0.10; % Drift ratio at which we are calling the model essentially collapsed
analysis.damping = 'rayleigh';
analysis.damp_ratio = 0.05;
analysis.joint_model = 1;
analysis.joint_explicit = 0;
analysis.additional_elements = 1;

% GM options
analysis.gm_set = 'FEMA_far_field';
analysis.run_sa_stripes = 1;
analysis.run_z_motion = 0;
analysis.scale_method = '2D'; % 'maxdir' or 'geomean' or '2D'
analysis.run_hazards = 0;
% analysis.rps2run = [43, 72, 108, 224, 475, 750, 975, 2475];
% analysis.rps2run = [475];
analysis.run_dbe = 1;
analysis.de_levels = [0.2 0.4 0.6 0.8 1 1.25 1.5];
analysis.run_mce = 0;
analysis.mce_levels = [1];
% analysis.scale_increment = 0.25;
% analysis.sa_stripes = [0.2 0.4];

% Nonlinear options
analysis.nonlinear_type = 'lumped'; % lumped or fiber
analysis.stories_nonlinear = inf;
analysis.stories_nonlinear_low = 0;
analysis.elastic_beams = 0;
analysis.fiber_walls = 0;
analysis.hinge_stiff_mod = 10;

% Solution algorithm options
analysis.solution_algorithm = 1;
analysis.initial_timestep_factor = 1;
analysis.algorithm = 'KrylovNewton';
analysis.integrator = 'Newmark 0.5 0.25';

% Run options
if ~isfield(analysis,'clear_existing_data')
    analysis.clear_existing_data = 1;
end
analysis.model_build_remote = 0;
analysis.summit = 0;
analysis.run_parallel = 1;
analysis.opensees_SP = 0;
analysis.run_ida = 1;
analysis.post_process_ida = 1;
analysis.create_fragilities = 0;
analysis.general_ida = 0;

% Visualization and Output options
% analysis.plot_ida = 0;
% analysis.detialed_post_process = 0;
analysis.write_xml = 1;
analysis.play_movie = 0;
analysis.movie_scale = 10;
analysis.suppress_outputs = 1;
if ~isfield(analysis,'simple_recorders')
    analysis.simple_recorders = 1;
end
analysis.hinge_stories_2_plot = 0;

end

