clear all
close all
clc

import usgs.*

% cmap1 = [0, 94, 162
%          37, 122, 45]/255;
     
% LA Site
site.lat = 34.05;
site.lng = -118.25;
site.vs30 = 537;
rps2run = [100 475 975 2475 4975];
[sa_spectra, sa_periods] = fn_call_USGS_hazard_API('E2014', site.lat, site.lng, site.vs30, 1./rps2run);
sa_1_LA = interp1(sa_periods,sa_spectra,1);

% SEA Site
site.lat = 47.6;
site.lng = -122.3;
site.vs30 = 537;
rps2run = [100 475 975 2475 4975];
[sa_spectra, sa_periods] = fn_call_USGS_hazard_API('E2014', site.lat, site.lng, site.vs30, 1./rps2run);
sa_1_SEA = interp1(sa_periods,sa_spectra,1);


% Write data as txt file
tab_dir = 'C:\Users\dtc2\OneDrive - NIST\NIST\Dissertation Archetype Study\Writeup\Tables';
fid = fopen([tab_dir filesep 'gm_table.tex'], 'w');
for i = 1:length(rps2run)
    if i == length(rps2run)
        fprintf(fid,'%i & %0.2f & %0.2f', rps2run(i), sa_1_LA(i), sa_1_SEA(i));
    else
        fprintf(fid,'%i & %0.2f & %0.2f \\\\\n', rps2run(i), sa_1_LA(i), sa_1_SEA(i));
    end
end
fclose(fid);



% plt = loglog(sa_periods,sa_spectra,'color',cmap1(1,:));
% plt(1).DisplayName = '100-year';
% plt(2).DisplayName = '475-year';
% plt(3).DisplayName = '975-year';
% plt(4).DisplayName = '2475-year';
% plt(5).DisplayName = '4975-year';
% xlabel('T (s)')
% ylabel('Sa (g)')
% xlim([0.1 3])
% ylim([0.01, 5])
% box on
% legend('location','southwest')


% 
% analysis.gm_set = 'FEMA_far_field';
% gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);
% 
% period = 1;
% targ_sa = sa_spectra(sa_periods == period);
% 
% % Run for each scale factor
% for i = 1:height(gm_set_table)
%     % Define Ground Motion Data
%     gm = gm_set_table(i,:);
%     eq_dir = {['ground_motions' '/' analysis.gm_set '/' gm.eq_name{1}]};
%     eq_name = {[gm.eq_name{1} '.tcl']};
% 
%     % Load spectral info and scale ground motion as geomean or pair
%     spectra_table = readtable([eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
%     sa_gm = interp1(spectra_table.period,spectra_table.psa_5,period);
%     scale_factor = targ_sa / sa_gm;
%     
%     % Save gm spectra
%     gm_spectra{i}.sa = spectra_table.psa_5 * scale_factor;
%     gm_spectra{i}.T = spectra_table.period;
% end
% 
% % Plot 
% for i = 1:height(gm_set_table)
%     semilogy(gm_spectra{i}.T,gm_spectra{i}.sa,'color',[0.5 0.5 0.5],'linewidth',0.5)
%     hold on
% end
% plot(sa_periods,sa_spectra,'k','linewidth',1.5)
% 
% xlabel('T (s)')
% ylabel('Sa (g)')
% xlim([0.1 3])
% box on
% grid on
% % ax = gca;
% % ax.YScale = 'log';