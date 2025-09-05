function [ eigen ] = main_eigen_analysis( analysis, read_dir_model )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import Packages
import opensees.*
import opensees.write_tcl.*
import file_exchange.*

% Create Write Directory
% TCL or Opensees does not like filesep command on windows, therefore must manually define forward slash seperators
write_dir_opensees = [strrep(analysis.out_dir,'\','/') '/eigen_analysis']; 
fn_make_directory( write_dir_opensees )

% Load Model Data
model = readtable([read_dir_model filesep 'model.csv'],'ReadVariableNames',true);
node = readtable([read_dir_model filesep 'node.csv'],'ReadVariableNames',true);
element = readtable([read_dir_model filesep 'element.csv'],'ReadVariableNames',true);
ele_props_table = readtable([read_dir_model filesep 'element_table.csv'],'readVariableNames',true);
story = readtable([read_dir_model filesep 'story.csv'],'ReadVariableNames',true);
joint = readtable([read_dir_model filesep 'joint.csv'],'ReadVariableNames',true);
hinge = readtable([read_dir_model filesep 'hinge.csv'],'ReadVariableNames',true);

%% Run Opensees for Eigen Analysis
% Write TCL files
analysis.horizontal_mass_only = 1;
[ ~ ] = fn_define_model( write_dir_opensees, node, element, joint, hinge, analysis, model.dimension, story, [], model, ele_props_table );
primary_nodes = node.id(node.primary_story == 1 & node.story > 0);
[ num_modes ] = fn_eigen_analysis( write_dir_opensees, primary_nodes', max(story.id), analysis, model.dimension);
fn_define_eignen_run_script( write_dir_opensees )

% Run Opensees
fprintf('Running Opensees Eigen Analysis\n')
main_run_opensees( write_dir_opensees, analysis )

% Postprocess OS data
periods = dlmread([write_dir_opensees filesep 'period.txt']);
model.T1_x = periods(1);

% Grab mode shapes
% hold on
mode_shapes = [];
for i = 1:num_modes
    xml_read = xmlread([write_dir_opensees filesep 'mode_shape_' num2str(i) '.xml']);
    [ xml_struct ] = xml2struct( xml_read );
    mode_raw = str2num(xml_struct.OpenSees.Data.Text);
%     mode_raw = dlmread([write_dir_opensees filesep 'mode_shape_' num2str(i) '.txt']);

    idx = 1:3:model.num_stories*3;
    mode = mode_raw(idx)';
    mode_shapes(:,i) = mode / mode(end);
%     plot(mode_norm(:,i),1:model.num_stories)
end

%% Save results
% periods
eigen.periods = periods';
writetable(model,[read_dir_model filesep 'model.csv'])

% mode shapes
eigen.mode_shapes = mode_shapes;
save([write_dir_opensees filesep 'mode_shapes.mat'],'mode_shapes')


end

