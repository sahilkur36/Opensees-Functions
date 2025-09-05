function [ ] = fn_pull_EDPs( main_dir, gm_set_table, story, analysis )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

disp('Collecting EDP data for FEMA P-58 Assessment ...')
id = 0;
idh = 0;
idr = [];
pfa = [];
col_rot = table;
bm_rot = table;
wall_rot = table;
for gm = 1:height(gm_set_table)
    gm_set_id = gm_set_table.set_id(gm);
    gm_pair_id = gm_set_table.pair(gm);
    alt_gm_pair_id = gm_set_table.pair(gm_set_table.set_id == gm_set_id & gm_set_table.pair ~= gm_pair_id);
    gm_dir     = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_id) '_' num2str(gm_pair_id)];
    alt_gm_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_id) '_' num2str(alt_gm_pair_id)];
    sa_folders = dir([gm_dir filesep 'Sa_*']);
    for s = 1:length(sa_folders)
        % Load data
        outputs_dir = [gm_dir '/' sa_folders(s).name];
        outputs_dir_alt_dir = [alt_gm_dir '/' sa_folders(s).name];
        outputs_file = [outputs_dir filesep 'summary_results.mat'];
        outputs_file_alt_dir = [outputs_dir_alt_dir filesep 'summary_results.mat'];
        story_file = [outputs_dir filesep 'story_analysis.mat'];
        hinge_file = [outputs_dir filesep 'hinge_analysis.mat'];
        ele_file = [outputs_dir filesep 'element_analysis.mat'];
        sa_val = round(str2double(strrep(strrep(sa_folders(s).name,'Sa_',''),'_','.')),3);
%             rp_val = rps2run(min(abs(analysis.sa_stripes - sa_val)) == abs(analysis.sa_stripes - sa_val));

        % Load Data
        if exist(outputs_file,'file')
            load(outputs_file)
            alt_dir_out = load(outputs_file_alt_dir);
        else
            error('NEED TO EXPLICITLY HANDLE MISSING SUMMARY TABLE')
        end

        % Define collapse based on two direction of motion
        building_collapse = summary.collapse | alt_dir_out.summary.collapse;
%             building_collapse = summary.collapse;

        % Component response
        idh = idh + 1;
        for n = 1:height(story)
            col_rot.sa(idh,1) = sa_val;
            col_rot.direction(idh,1) = gm_set_table.pair(gm);
            col_rot.gm(idh,1) = gm_set_table.set_id(gm);
            bm_rot.sa(idh,1) = sa_val;
            bm_rot.direction(idh,1) = gm_set_table.pair(gm);
            bm_rot.gm(idh,1) = gm_set_table.set_id(gm);
            wall_rot.sa(idh,1) = sa_val;
            wall_rot.direction(idh,1) = gm_set_table.pair(gm);
            wall_rot.gm(idh,1) = gm_set_table.set_id(gm);

            if analysis.nonlinear && exist(hinge_file,'file') && ~building_collapse
                load(hinge_file)
                max_col_deform = max(hinge.tot_deform(strcmp(hinge.ele_type,'column') & hinge.story == n));
                if isempty(max_col_deform)
                    col_rot.(['story_' num2str(n)])(idh,1) = NaN;
                else
                    col_rot.(['story_' num2str(n)])(idh,1) = max_col_deform;
                end
                max_beam_deform = max(hinge.tot_deform(strcmp(hinge.ele_type,'beam') & hinge.story == n));
                if isempty(max_beam_deform)
                    bm_rot.(['story_' num2str(n)])(idh,1) = NaN;
                else
                    bm_rot.(['story_' num2str(n)])(idh,1) = max_beam_deform;
                end
                max_wall_deform = max(hinge.tot_deform(strcmp(hinge.ele_type,'wall') & hinge.story == n));
                if isempty(max_wall_deform)
                    wall_rot.(['story_' num2str(n)])(idh,1) = NaN;
                else
                    wall_rot.(['story_' num2str(n)])(idh,1) = max_wall_deform;
                end
            elseif ~analysis.nonlinear && exist(ele_file,'file') && ~building_collapse
                ele_results = load(ele_file);
                max_col_deform = max(ele_results.element.max_rot(strcmp(ele_results.element.type,'column') & ~ele_results.element.rigid & ele_results.element.story == n));
                if isempty(max_col_deform)
                    col_rot.(['story_' num2str(n)])(idh,1) = NaN;
                else
                    col_rot.(['story_' num2str(n)])(idh,1) = max_col_deform;
                end
                max_beam_deform = max(ele_results.element.max_rot(strcmp(ele_results.element.type,'beam') & ~ele_results.element.rigid & ele_results.element.story == n));
                if isempty(max_beam_deform)
                    wall_rot.(['story_' num2str(n)])(idh,1) = NaN;
                else
                    bm_rot.(['story_' num2str(n)])(idh,1) = max_beam_deform;
                end
                max_wall_deform = max(ele_results.element.max_rot(strcmp(ele_results.element.type,'wall') & ~ele_results.element.rigid & ele_results.element.story == n));
                if isempty(max_wall_deform)
                    wall_rot.(['story_' num2str(n)])(idh,1) = NaN;
                else
                    wall_rot.(['story_' num2str(n)])(idh,1) = max_wall_deform;
                end
            else
                col_rot.(['story_' num2str(n)])(idh,1) = NaN;
                bm_rot.(['story_' num2str(n)])(idh,1) = NaN;
                wall_rot.(['story_' num2str(n)])(idh,1) = NaN;
            end
        end

        % EDPs
        if building_collapse == 0
            if exist(story_file,'file')
                load(story_file)

%                     for d = 1:2
                    id = id + 1;
                    % Formulate the IDR table for the P-58 assessment - X direction
                    idr.sa(id,1) = sa_val;
%                         idr.rp(id,1) = rp_val;
                    idr.direction(id,1) = gm_set_table.pair(gm);
%                         idr.direction(id,1) = d;
                    idr.gm(id,1) = gm_set_table.set_id(gm);
%                         idr.gm(id,1) = gm_set_table.id(gm);
                    for n = 1:height(story)
                        if building_collapse % set collapse values = NaN
                            idr.(['story_' num2str(n)])(id,1) = NaN;
                        else
                            idr.(['story_' num2str(n)])(id,1) = story.max_drift_x(n);
                        end
                    end
%                     end

                    % Formulate the PFA table for the P-58 assessment - X direction
                    pfa.sa(id,1) = sa_val;
%                         pfa.rp(id,1) = rp_val;
%                     for d = 1:2
                    pfa.direction(id,1) = gm_set_table.pair(gm);
%                         pfa.direction(id,1) = d;
                    pfa.gm(id,1) = gm_set_table.set_id(gm);
%                         pfa.gm(id,1) = gm_set_table.id(gm);
                    pfa.floor_1(id,1) = summary.pga_x;
                    for n = 1:height(story)
                        if building_collapse % set collapse values = NaN
                            pfa.(['floor_' num2str(n+1)])(id,1) = NaN;
                        else
                            pfa.(['floor_' num2str(n+1)])(id,1) = story.max_accel_x(n);
                        end
                    end
%                     end

%                     % Pull z direction response - need to update this if
%                     I actually run
%                     if analysis.run_z_motion
%                         id = id + 1;
%                         % Formulate the IDR table for the P-58 assessment - Z direction 
%                         idr.sa(id,1) = sa_val;
%                         idr.direction(id,1) = 2;
%                         idr.gm(id,1) = gm;%gm_set_table.eq_name{gm};
%                         for n = 1:height(story)
%                             if summary.collapse > 0 % set collapse values = NaN
%                                 idr.(['story_' num2str(n)])(id,1) = NaN;
%                             else
%                                 idr.(['story_' num2str(n)])(id,1) = story.max_drift_z(n);
%                             end
%                         end
% 
%                         % Formulate the PFA table for the P-58 assessment - Z direction 
%                         pfa.sa(id,1) = sa_val;
%                         pfa.direction(id,1) = 2;
%                         pfa.gm(id,1) = gm;%gm_set_table.eq_name{gm};
%                         pfa.floor_1(id,1) = summary.pga_z;
%                         for n = 1:height(story)
%                             if summary.collapse > 0 % set collapse values = NaN
%                                 pfa.(['floor_' num2str(n+1)])(id,1) = NaN;
%                             else
%                                 pfa.(['floor_' num2str(n+1)])(id,1) = story.max_accel_z(n);
%                             end
%                         end
%                     end
            else
                error('NEED TO EXPLICITLY HANDLE MISSING STORY TABLE')
            end
        end
    end
end

% Convert to tables and sort
idr_table = struct2table(idr);
pfa_table = struct2table(pfa);
idr_table_sort = [];
pfa_table_sort = [];
sa_vals = unique(idr_table.sa);
for i = 1:length(sa_vals)
    for d = 1:2
        filt = idr_table.sa == sa_vals(i) & idr_table.direction == d;
        idr_table_sect = sortrows(idr_table(filt,:),3);
        idr_table_sort = [idr_table_sort; idr_table_sect];
        pfa_table_sect = sortrows(pfa_table(filt,:),3);
        pfa_table_sort = [pfa_table_sort; pfa_table_sect];
    end
end

% Write tables to save directory
write_dir = [main_dir '/' 'IDA' '/' 'EDPs'];
mkdir(write_dir)
writetable(idr_table_sort,[write_dir filesep 'idr.csv'])
writetable(pfa_table_sort,[write_dir filesep 'pfa.csv'])
writetable(col_rot,[write_dir filesep 'col_deform.csv'])
writetable(bm_rot,[write_dir filesep 'beam_deform.csv'])
writetable(wall_rot,[write_dir filesep 'wall_deform.csv'])

end

