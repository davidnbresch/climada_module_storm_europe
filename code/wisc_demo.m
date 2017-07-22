% wisc_demo
% climada template
% MODULE:
%   storm_europe
% NAME:
%   wisc_demo
% PURPOSE:
%   batch job to run WISC demo (key function called: wisc_hazard_set)
%
%   FIRST, set PARAMETERS in code below, such as the folder with the WISC
%   netCDF storm footprints
%
%   see comments in code below. For speedup in subsequent calls, the WISC
%   hazard event sets are stored.
%
% CALLING SEQUENCE:
%   wisc_demo
% EXAMPLE:
%   wisc_demo
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to stdout and figures stores as .png files
%   the WISC hazard sets are stored into the climada hazards folder
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20170721, initial
% David N. Bresch, david.bresch@gmail.com, 20170722, Octave version (few plots, for speedup)
% David N. Bresch, david.bresch@gmail.com, 20170722, EM-DAT comparison added
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
% this section defines all parameters
%
% define the TEST country
country_name='GBR'; % name like 'United Kingdom' or ISO3 code like 'GBR'
%country_name='FRA'; % name like 'United Kingdom' or ISO3 code like 'GBR'
%
% local folder with the WISC netCDF storm footprints
wisc_dir='/Users/bresch/polybox/WISC';
%
% local folder to write the figures
fig_dir ='/Users/bresch/Desktop/WISC';
%
CAGR=0.02; % compound annual growth rate to inflate EM-DAT damages
%
% whether we plot all hazard statistics
plot_hazard=0; % default=1

% obtain country name and ISO3 code
[country_name,country_ISO3]=climada_country_name(country_name);
country_ISO3_title=[country_name ' (' country_ISO3 ')'];
ISO3_country_name=[country_ISO3 '_' strrep(country_name,' ','')];


% create the WISC hazard event sets (if it does not exist)
% ---------------------------------

hazard_era20c_file=[climada_global.hazards_dir filesep 'WISC_era20c_eur_WS.mat'];
if ~exist(hazard_era20c_file,'file')
    % generate the hazard event set from netCDF footprints
    % ====================================================
    hazard_era20c=wisc_hazard_set([wisc_dir filesep 'fp_era20c_*.nc'],0,'WISC_era20c_eur_WS');
else
    % load the previously generated hazard event set
    hazard_era20c=climada_hazard_load(hazard_era20c_file);
end

hazard_eraint_file=[climada_global.hazards_dir filesep 'WISC_eraint_eur_WS.mat'];
if ~exist(hazard_eraint_file,'file')
    hazard_eraint=wisc_hazard_set([wisc_dir filesep 'fp_eraint_*.nc'],0,'WISC_eraint_eur_WS');
else
    hazard_eraint=climada_hazard_load(hazard_eraint_file);
end


% show some hazard statistics
% ---------------------------

if ~climada_global.octave_mode && plot_hazard
    figure;res=climada_hazard_plot(hazard_era20c,0); % max intensity at each centroid
    title(['ERA_{20c} ' res.title_str]) % just to show also source of footprints
    saveas(gcf,[fig_dir filesep 'ERA20c_max_intens.png'],'png');
    figure;res=climada_hazard_plot(hazard_era20c,-1); % strongest storm
    title(['ERA_{20c} ' res.title_str]) % just to show also source of footprints
    saveas(gcf,[fig_dir filesep 'ERA20c_max_storm.png'],'png');
    figure;res=climada_hazard_plot(hazard_eraint,0); % max intensity at each centroid
    title(['ERA_{int} ' res.title_str]) % just to show also source of footprints
    saveas(gcf,[fig_dir filesep 'ERAint_max_intens.png'],'png');
    figure;res=climada_hazard_plot(hazard_eraint,-1); % strongest storm
    title(['ERA_{int} ' res.title_str]) % just to show also source of footprints
    saveas(gcf,[fig_dir filesep 'ERAint_max_storm.png'],'png');
else
    fprintf('Octave: plots suppressed for speedup\n');
end

% create the (default, 10x10km) asset base
% ----------------------------------------

entity_file=[climada_global.entities_dir filesep 'WISC_' country_ISO3 '_' strrep(country_name,' ','') '_10x10km.mat'];
if ~exist(entity_file,'file')
    % generate the asset data for requested country on 10x10km resolution
    % ===================================================================
    % (includes the default damage function for storm Europe, too)
    entity=climada_entity_country(country_ISO3);
    entity=climada_assets_encode(entity,hazard_era20c);
    save(entity_file,'entity','-v7.3'); % save as hdf5
else
    % load previously generated assets
    entity=climada_entity_load(entity_file);
end
figure;climada_entity_plot(entity); % plot it
title(['Assets for ' country_ISO3_title ' (10x10km)'])
saveas(gcf,[fig_dir filesep ISO3_country_name '_assets.png'],'png');


% calculate damages for all events for asset base
% ===============================================
% EDS = event damage set
clear EDS % just in case we call this batch code again
EDS(1)=climada_EDS_calc(entity,hazard_era20c);
EDS(2)=climada_EDS_calc(entity,hazard_eraint); % assume same hazard resolution


% plot the exceedance damage frequency curves (DFC)
% -------------------------------------------------

figure;[DFC,fig,legend_str]=climada_EDS_DFC(EDS);title(['Assets for ' country_ISO3_title])
saveas(gcf,[fig_dir filesep ISO3_country_name '_DFC.png'],'png');


% add EM-DAT
em_data=emdat_read('',country_name,'WS',1,1,CAGR); % 2nd last parameter =1 for verbose
if ~isempty(em_data)
    hold on
    % next would be attemps to automatically adjust:
    % [adj_EDS1,climada2emdat_factor_weighted] = cr_EDS_emdat_adjust(EDS(1));
    % if abs(climada2emdat_factor_weighted-1)>10*eps
    %     adj_DFC1=climada_EDS2DFC(adj_EDS1);
    %     plot(adj_DFC1.return_period,adj_DFC1.damage,':y','LineWidth',1);
    %     legend_str{end+1}='era20c EM-DAT adjusted';
    % end
    % [adj_EDS2,climada2emdat_factor_weighted] = cr_EDS_emdat_adjust(EDS(2));
    % if abs(climada2emdat_factor_weighted-1)>10*eps
    %     adj_DFC2=climada_EDS2DFC(adj_EDS2);
    %     plot(adj_DFC2.return_period,adj_DFC2.damage,':r','LineWidth',1);
    %     legend_str{end+1}='eraint EM-DAT adjusted';
    % end
    
    plot(em_data.DFC.return_period,em_data.DFC.damage,'dg');
    legend_str{end+1} = em_data.DFC.annotation_name;
    plot(em_data.DFC.return_period,em_data.DFC_orig.damage,'og');
    legend_str{end+1} = em_data.DFC_orig.annotation_name;
    legend(legend_str,'Location','NorthWest'); % show legend
    saveas(gcf,[fig_dir filesep ISO3_country_name '_DFC_emdat.png'],'png');
end % em_data


% show the most damaging storms
% -----------------------------

if ~climada_global.octave_mode
    [~,max_pos]=max(EDS(1).damage);
    figure;res=climada_hazard_plot(hazard_era20c,max_pos); % max intensity at each centroid
    title(['ERA_{20c} ' res.yyyymmdd_str ' most damaging for ' country_ISO3_title]);
    xlabel(sprintf('damage %2.0f mio USD',EDS(1).damage(max_pos)/1e6));ylabel('');
    saveas(gcf,[fig_dir filesep ISO3_country_name '_ERA20c_max_damage.png'],'png');
    
    [~,max_pos]=max(EDS(2).damage);
    figure;res=climada_hazard_plot(hazard_eraint,max_pos); % max intensity at each centroid
    title(['ERA_{int} ' res.yyyymmdd_str ' most damaging for ' country_ISO3_title]);
    xlabel(sprintf('damage %2.0f mio USD',EDS(2).damage(max_pos)/1e6));ylabel('');
    saveas(gcf,[fig_dir filesep ISO3_country_name '_ERAint_max_damage.png'],'png');
else
    fprintf('Octave: plots suppressed for speedup\n');
end