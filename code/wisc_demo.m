% wisc_demo
% wisc_demo
% MODULE:
%   storm_europe
% NAME:
%   wisc_demo
% PURPOSE:
%   batch job to run WISC demo (key function called: wisc_hazard_set)
%
%   FIRST, set PARAMETERS in code below, such as the folder with the WISC
%   netCDF storm footprints. Visit https://wisc.climate.copernicus.eu/wisc/#/help/products 
%   and download C3S_WISC_FOOTPRINT_NETCDF_0100.tgz, then edit the wisc_dir below 
%
%   see comments in code below. For speedup in subsequent calls, the WISC
%   hazard event sets are stored. If run in Ocatve, only damage frequency
%   curve (DFC) is shown, as plotting of hazard and entity takes a lot of
%   time (Octave is much slower in plotting than MATLAB).
%
%   SPECIAL: if there is a WS hazard set for the specified country, the
%   code does also encode to and compare with.
%   If you further provide a hazard set as hazard_cmp (being loaded before
%   you call wisc_demo), the code does also encode to and compare with.
%
%   See show_BOTH_separate in PARAMETERS section in code to run both parts
%   (era20 and eraint) separately
%
%   check:
%   fp_era20c_1972111300_177_0
%   fp_era20c_1972111309_177_0
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
% David N. Bresch, david.bresch@gmail.com, 20170725, storm module comparison added
% David N. Bresch, david.bresch@gmail.com, 20170726, combined view added
% David N. Bresch, david.bresch@gmail.com, 20170730, Octave use improved
% David N. Bresch, david.bresch@gmail.com, 20170731, Octave printing
% David N. Bresch, david.bresch@gmail.com, 20171004, hint to WISC data source, paths relative
% David N. Bresch, david.bresch@gmail.com, 20171108, combined full hazard event set, last (duplicate) era20 event excluded
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
% this section defines all parameters
%
% very special case to show the two sub-sets
show_BOTH_separate=0; % default=0
%
% define the TEST country or -ies (at least one)
%country_names={'GBR','IRL','DEU','FRA','DNK','NLD','BEL','CHE','AUT','ESP'}; % name like 'United Kingdom' or ISO3 code like 'GBR'
country_names={'GBR','DEU','FRA'}; % name like 'United Kingdom' or ISO3 code like 'GBR'
%country_names={'GBR'}; % name like 'United Kingdom' or ISO3 code like 'GBR'
%country_names={'GBR','IRL','DEU','FRA','DNK','NLD','BEL','CHE','ESP'}; % name like 'United Kingdom' or ISO3 code like 'GBR'
%country_names={'GBR','DEU','FRA'}; % name like 'United Kingdom' or ISO3 code like 'GBR'
%country_names={'GBR'}; % name like 'United Kingdom' or ISO3 code like 'GBR'
%country_names={'Germany','Switzerland'}; % name like 'United Kingdom' or ISO3 code like 'GBR'
%
% local folder with the WISC netCDF storm footprints
%wisc_dir='/Users/bresch/polybox/WISC/footprints';
%wisc_dir='C:\shortpaths\Documents\WISC\polybox_clone\footprints';
%wisc_dir='D:\Documents_DATA\WISC_data_20170918\Storm Footprints';
wisc_dir=[climada_global.data_dir filesep 'WISC' filesep 'C3S_WISC_FOOTPRINT_NETCDF_0100'];
if ~isdir(wisc_dir)
    fprintf('Please locate the WISC netCDF files folder:\n');
    wisc_dir = uigetdir(climada_global.data_dir, 'Locate WISC netCDF files folder:');
    if isequal(filename,0) || isequal(pathname,0)
        fprintf('No WISC netCDF file folder selected\n');
        fprintf(['> visit <a href="https://wisc.climate.copernicus.eu/wisc/#/help/products">'...
            'https://wisc.climate.copernicus.eu/wisc/#/help/products</a>\n',...
            '  and download C3S_WISC_FOOTPRINT_NETCDF_0100.tgz\n']);
        return
    end
end % ~isdir(wisc_dir)
%
% local folder to write the figures
fig_dir =[climada_global.results_dir filesep 'WISC'];
if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
fig_ext ='png';
%
% whether we plot all hazard statistics
plot_hazard=0; % default=1
plot_entity=0; % default=1
%
% whether we compare with the default climada (hazard,entity,vulnerability)
standard_storm_module=0; % default=0
%
climada_global_save_file_version=climada_global.save_file_version; % store
climada_global.save_file_version='-v7'; % NOT -v7.3 to allow reading in Octave

if climada_global.octave_mode
    % since Octave is VERY slow in plotting, skip this unless requested
    plot_hazard=1;
    plot_entity=1;
    sok=1; % set this =0 to skip saving to disk (if troubles)
    fig_ext ='eps'; % saving to .png is not possible
else
    sok=1; % saving to .png is ok
end

% create the WISC hazard event sets (if not existing)
% ---------------------------------

hazard_wisc_file=[climada_global.hazards_dir filesep 'WISC_eur_WS.mat'];
if exist('hazard_wisc','var')
    fprintf('NOTE: using previously loaded hazard_wisc\n');
elseif exist(hazard_wisc_file,'file')
    fprintf('loading %s\n',hazard_wisc_file);
    hazard_wisc=climada_hazard_load(hazard_wisc_file);
else
    if exist('hazard_era20c','var') && exist('hazard_eraint','var')
        fprintf('NOTE: using previously loaded hazard_era20c and hazard_eraint\n');
    else
        
        hazard_era20c_file=[climada_global.hazards_dir filesep 'WISC_era20c_eur_WS.mat'];
        if ~exist(hazard_era20c_file,'file')
            % generate the hazard event set from netCDF footprints
            % ====================================================
            hazard_era20c=wisc_hazard_set([wisc_dir filesep 'fp_era20c_*.nc'],0,'WISC_era20c_eur_WS');
        else
            fprintf('NOTE: loading previously calculated hazard_era20c ...\n');
            hazard_era20c=climada_hazard_load(hazard_era20c_file,1);
        end
        
        hazard_eraint_file=[climada_global.hazards_dir filesep 'WISC_eraint_eur_WS.mat'];
        if ~exist(hazard_eraint_file,'file')
            hazard_eraint=wisc_hazard_set([wisc_dir filesep 'fp_eraint_*.nc'],0,'WISC_eraint_eur_WS');
        else
            fprintf('NOTE: loading previously calculated hazard_eraint ...\n');
            hazard_eraint=climada_hazard_load(hazard_eraint_file,1);
        end
        
        fprintf('NOTE: combining hazard_era20c and hazard_eraint into one hazard set\n');

        hazard=climada_hazard_merge(hazard_era20c,hazard_eraint,'events');
        hazard=rmfield(hazard,'lonlat_size');
        save(hazard_wisc_file,'hazard','-v7.3');
        hazard_wisc=hazard;clear hazard_wisc
        
    end
end

% show some hazard statistics
% ---------------------------

if plot_hazard
    % one could also use climada_hazard_plot, which contours the plot, but slower
    if show_BOTH_separate
        figure;res=climada_hazard_plot_nogrid(hazard_era20c, 0,1); % max intensity at each centroid
        title(['ERA 20c ' res.title_str]) % just to show also source of footprints
        if sok,saveas(gcf,[fig_dir filesep 'ERA20c_max_intens.' fig_ext],fig_ext);end
        figure;res=climada_hazard_plot_nogrid(hazard_era20c,-1,1); % strongest storm
        title(['ERA 20c ' res.title_str]) % just to show also source of footprints
        if sok,saveas(gcf,[fig_dir filesep 'ERA20c_max_storm.' fig_ext],fig_ext);end
        figure;res=climada_hazard_plot_nogrid(hazard_eraint, 0,1); % max intensity at each centroid
        title(['ERA int ' res.title_str]) % just to show also source of footprints
        if sok,saveas(gcf,[fig_dir filesep 'ERAint_max_intens.' fig_ext],fig_ext);end
        figure;res=climada_hazard_plot_nogrid(hazard_eraint,-1,1); % strongest storm
        title(['ERA int ' res.title_str]) % just to show also source of footprints
        if sok,saveas(gcf,[fig_dir filesep 'ERAint_max_storm.' fig_ext],fig_ext);end
    else
        figure;res=climada_hazard_plot_nogrid(hazard_wisc, 0,1); % max intensity at each centroid
        title(['WISC ' res.title_str]) % just to show also source of footprints
        if sok,saveas(gcf,[fig_dir filesep 'WISC_max_intens.' fig_ext],fig_ext);end
        figure;res=climada_hazard_plot_nogrid(hazard_wisc,-1,1); % strongest storm
        title(['WISC ' res.title_str]) % just to show also source of footprints
        if sok,saveas(gcf,[fig_dir filesep 'WISC_max_storm.' fig_ext],fig_ext);end
    end
end

% loop over countries for a first impression

entity_combined=[];EDS_era20c=[];EDS_eraint=[];EDS_std_comb=[];EDS_wisc=[];  % init
em_data_ET_tot=[];em_data_ET_tot.year=[];em_data_ET_tot.damage=[];em_data_ET_tot.damage_orig=[]; % init
em_data_ST_tot=[];em_data_ST_tot.year=[];em_data_ST_tot.damage=[];em_data_ST_tot.damage_orig=[]; % init

for country_i=1:length(country_names)
    
    country_name=country_names{country_i};
    
    % obtain country name and ISO3 code
    [country_name,country_ISO3]=climada_country_name(country_name);
    country_ISO3_title=[country_name ' (' country_ISO3 ')'];
    ISO3_country_name=[country_ISO3 '_' strrep(country_name,' ','')];
    
    fprintf('\n*** processing %s: ***\n',country_ISO3_title);
    
    % create the (default, 10x10km) asset base
    % ----------------------------------------
    
    entity_file=[climada_global.entities_dir filesep 'WISC_' country_ISO3 '_' strrep(country_name,' ','') '_10x10km.mat'];
    if ~exist(entity_file,'file')
        % generate the asset data for requested country on 10x10km resolution
        % ===================================================================
        % (includes the default damage function for storm Europe, too)
        entity=climada_entity_country(country_ISO3);
        entity=climada_assets_encode(entity,hazard_era20c);
        save(entity_file,'entity',climada_global.save_file_version);
    else
        % load previously generated assets
        entity=climada_entity_load(entity_file,1);
        save(entity_file,'entity','-v7'); % for Octave compatibility

    end
    if plot_entity
        figure;climada_entity_plot(entity); % plot it
        title(['Assets for ' country_ISO3_title ' (10x10km)'])
        if sok,saveas(gcf,[fig_dir filesep ISO3_country_name '_assets.' fig_ext],fig_ext);end
    end % plot_entity
    
    % add to combined entity
    entity_combined=climada_entity_combine(entity_combined,entity);
    
    
    % calculate damages for all events for asset base
    % ===============================================
    % EDS = event damage set
    clear EDS % just in case we call this batch code again
    if show_BOTH_separate
        EDS(1)=climada_EDS_calc(entity,hazard_era20c);
        EDS(2)=climada_EDS_calc(entity,hazard_eraint); % assume same hazard resolution
        EDS(1).annotation_name='ERA 20c'; % as the name got long
        EDS(2).annotation_name='ERA int';
    else
        EDS(1)=climada_EDS_calc(entity,hazard_wisc);
        EDS(1).annotation_name='WISC'; % as the name got long
    end
    
    
    % plot the exceedance damage frequency curves (DFC)
    % -------------------------------------------------
    
    figure;[~,~,legend_str,legend_handle]=climada_EDS_DFC(EDS);title(['Assets for ' country_ISO3_title])
    if sok,saveas(gcf,[fig_dir filesep ISO3_country_name '_DFC.' fig_ext],fig_ext);end
    
    if ~climada_global.octave_mode
        % add EM-DAT (historic event data for comparison)
        % -----------------------------------------------
        % first all exratropical storms
        em_data_ET=emdat_read('',country_ISO3,'WS',1,1); % last parameter =1 for verbose
        [legend_str,legend_handle]=emdat_barplot(em_data_ET,'dm','om','EM-DAT extratrop',legend_str,legend_handle,'SouthEast');
        % second all storms (since classification ambiguous)
        em_data_ST=emdat_read('',country_ISO3,'-WS',1,1); % last parameter =1 for verbose
        [legend_str,legend_handle]=emdat_barplot(em_data_ST,'db','ob','EM-DAT storms',legend_str,legend_handle,'SouthEast');legend('boxoff')
        if sok,saveas(gcf,[fig_dir filesep ISO3_country_name '_DFC_emdat.' fig_ext],fig_ext);end
        
        % some data collection (for annual aggregate further below)
        em_data_ET_tot.year        = [em_data_ET_tot.year   ; em_data_ET.year];
        em_data_ET_tot.damage      = [em_data_ET_tot.damage ; em_data_ET.damage];
        em_data_ET_tot.damage_orig = [em_data_ET_tot.damage_orig ; em_data_ET.damage_orig];
        em_data_ET_tot.last_year   = em_data_ET.last_year;
        em_data_ST_tot.year        = [em_data_ST_tot.year   ; em_data_ST.year];
        em_data_ST_tot.damage      = [em_data_ST_tot.damage ; em_data_ST.damage];
        em_data_ST_tot.damage_orig = [em_data_ST_tot.damage_orig ; em_data_ST.damage_orig];
        em_data_ST_tot.last_year   = em_data_ST.last_year;
    end % climada_global.octave_mode
    
    % some data collection (for annual aggregate further below)
    if show_BOTH_separate
        EDS_era20c                 = climada_EDS_combine(EDS_era20c,EDS(1));
        EDS_eraint                 = climada_EDS_combine(EDS_eraint,EDS(2));
    else
        EDS_wisc                   = climada_EDS_combine(EDS_wisc,EDS(1));
    end
    
    % load the standard hazard set for this country (if exists)
    hazard_std_file=[climada_global.hazards_dir filesep ISO3_country_name '_eur_WS.mat'];
    if exist(hazard_std_file,'file') && standard_storm_module
        hazard_std=climada_hazard_load(hazard_std_file);
        entity_std=climada_assets_encode(entity,hazard_std);
        EDS_std=climada_EDS_calc(entity_std,hazard_std);
        hold on
        DFC_std=climada_EDS2DFC(EDS_std);
        RP100_pos=find(DFC_std.return_period==100);
        legend_handle(end+1) = plot(DFC_std.return_period(1:RP100_pos),DFC_std.damage(1:RP100_pos),'or');
        legend_str{end+1}    = 'standard storm module';
        legend(legend_handle,legend_str,'Location','SouthEast');legend('boxoff')
        xlim([0 100])
        if sok,saveas(gcf,[fig_dir filesep ISO3_country_name '_DFC_emdat_std.' fig_ext],fig_ext);end
        EDS_std_comb=climada_EDS_combine(EDS_std_comb,EDS_std);
    end % compare_with_storm_europe_module
    
    if exist('hazard_cmp','var') % see SPECIAL in header above
        entity_cmp=climada_assets_encode(entity,hazard_cmp);
        EDS_cmp=climada_EDS_calc(entity_cmp,hazard_cmp);
        hold on
        DFC_cmp=climada_EDS2DFC(EDS_cmp);
        legend_handle(end+1) = plot(DFC_cmp.return_period(1:RP100_pos),DFC_cmp.damage(1:RP100_pos),'xk');
        legend_str{end+1}    = 'comparison';
        legend(legend_handle,legend_str,'Location','SouthEast');
        if sok,saveas(gcf,[fig_dir filesep ISO3_country_name '_DFC_emdat_std_cmp.' fig_ext],fig_ext);end
    end
    
    % show the most damaging storms
    % -----------------------------
    
    if plot_hazard
        % one could also use climada_hazard_plot, which contours the plot, but slower
        if show_BOTH_separate
            [~,max_pos]=max(EDS(1).damage);
            figure;res=climada_hazard_plot_nogrid(hazard_era20c,max_pos,1); % max intensity at each centroid
            title(['ERA 20c ' res.yyyymmdd_str ' most damaging for ' country_ISO3_title]);
            xlabel(sprintf('damage %2.0f mio USD',EDS(1).damage(max_pos)/1e6));ylabel('');
            if sok,saveas(gcf,[fig_dir filesep ISO3_country_name '_ERA20c_max_damage.' fig_ext],fig_ext);end
            [~,max_pos]=max(EDS(2).damage);
            figure;res=climada_hazard_plot_nogrid(hazard_eraint,max_pos,1); % max intensity at each centroid
            title(['ERA int ' res.yyyymmdd_str ' most damaging for ' country_ISO3_title]);
            xlabel(sprintf('damage %2.0f mio USD',EDS(2).damage(max_pos)/1e6));ylabel('');
            if sok,saveas(gcf,[fig_dir filesep ISO3_country_name '_ERAint_max_damage.' fig_ext],fig_ext);end
        else
            [~,max_pos]=max(EDS(1).damage);
            figure;res=climada_hazard_plot_nogrid(hazard_wisc,max_pos,1); % max intensity at each centroid
            title(['WISC ' res.yyyymmdd_str ' most damaging for ' country_ISO3_title]);
            xlabel(sprintf('damage %2.0f mio USD',EDS(1).damage(max_pos)/1e6));ylabel('');
            if sok,saveas(gcf,[fig_dir filesep ISO3_country_name '_WISC_max_damage.' fig_ext],fig_ext);end
        end
    else
        fprintf('Note: hazard plots suppressed for speedup\n');
    end
    
end % country_i

if show_BOTH_separate
    EDS_era20c(1).annotation_name='ERA20c'; % as the name got long
    EDS_eraint(1).annotation_name='ERAint';
else
    EDS_wisc(1).annotation_name='WISC'; % as the name got long
end

% plot combined (pan-European) annual aggregate DFC
% -------------------------------------------------

if show_BOTH_separate
    YDS_era20c=climada_EDS2YDS(EDS_era20c,hazard_era20c);
    YDS_eraint=climada_EDS2YDS(EDS_eraint,hazard_eraint);
    figure;[~,~,legend_str,legend_handle]=climada_EDS_DFC(YDS_era20c,YDS_eraint);title('combined annual aggregate DFC')
    if sok,saveas(gcf,[fig_dir filesep 'combined_aggregate_DFC.' fig_ext],fig_ext);end
else
    YDS_wisc=climada_EDS2YDS(EDS_wisc,hazard_wisc);
    figure;[~,~,legend_str,legend_handle]=climada_EDS_DFC(YDS_wisc);title('combined annual aggregate DFC')
    if sok,saveas(gcf,[fig_dir filesep 'combined_aggregate_DFC.' fig_ext],fig_ext);end
end

% add EM-DAT (we need to construct the DFC here first)
% ----------------------------------------------------

% first for extratropical storms (ET)
if ~isempty(em_data_ET_tot.year)
    % create the annual aggregate (just summing up all events in each year)
    unique_years=unique(em_data_ET_tot.year);
    for year_i=1:length(unique_years)
        em_data_ET_tot.annual_damage(year_i)     =sum(em_data_ET_tot.damage(     em_data_ET_tot.year==unique_years(year_i)));
        em_data_ET_tot.annual_damage_orig(year_i)=sum(em_data_ET_tot.damage_orig(em_data_ET_tot.year==unique_years(year_i)));
    end
    
    % create the annual aggregate exceedence frequency curve
    em_data_ET_tot.frequency = ones(length(em_data_ET_tot.annual_damage),1)/(em_data_ET_tot.last_year-min(unique_years)+1);
    [em_data_ET_tot.DFC.damage,em_data_ET_tot.DFC.return_period]=...
        climada_damage_exceedence(em_data_ET_tot.annual_damage',     em_data_ET_tot.frequency',[],1);
    [em_data_ET_tot.DFC_orig.damage,em_data_ET_tot.DFC_orig.return_period]=...
        climada_damage_exceedence(em_data_ET_tot.annual_damage_orig',em_data_ET_tot.frequency',[],1);
    
    [legend_str,legend_handle]=emdat_barplot(em_data_ET_tot,'dm','om','EM-DAT extratrop',legend_str,legend_handle,'SouthEast');
    legend('boxoff')
    if sok,saveas(gcf,[fig_dir filesep 'combined_aggregate_DFC_emdat.' fig_ext],fig_ext);end
end % em_data

% and again for storms (ST)
if ~isempty(em_data_ST_tot.year)
    % create the annual aggregate (just summing up all events in each year)
    unique_years=unique(em_data_ST_tot.year);
    for year_i=1:length(unique_years)
        em_data_ST_tot.annual_damage(year_i)     =sum(em_data_ST_tot.damage(     em_data_ST_tot.year==unique_years(year_i)));
        em_data_ST_tot.annual_damage_orig(year_i)=sum(em_data_ST_tot.damage_orig(em_data_ST_tot.year==unique_years(year_i)));
    end
    
    % create the annual aggregate exceedence frequency curve
    em_data_ST_tot.frequency = ones(length(em_data_ST_tot.annual_damage),1)/(em_data_ST_tot.last_year-min(unique_years)+1);
    [em_data_ST_tot.DFC.damage,em_data_ST_tot.DFC.return_period]=...
        climada_damage_exceedence(em_data_ST_tot.annual_damage',     em_data_ST_tot.frequency',[],1);
    [em_data_ST_tot.DFC_orig.damage,em_data_ST_tot.DFC_orig.return_period]=...
        climada_damage_exceedence(em_data_ST_tot.annual_damage_orig',em_data_ST_tot.frequency',[],1);
    
    [legend_str,legend_handle]=emdat_barplot(em_data_ST_tot,'db','ob','EM-DAT storms',   legend_str,legend_handle,'SouthEast');
    legend('boxoff')
    if sok,saveas(gcf,[fig_dir filesep 'combined_aggregate_DFC_emdat.' fig_ext],fig_ext);end
end % em_data

if ~isempty(EDS_std_comb)
    hold on
    YDS_std_comb=climada_EDS2YDS(EDS_std_comb);
    YFC_std=climada_EDS2DFC(YDS_std_comb);
    RP100_pos=find(YFC_std.return_period==100);
    legend_handle(end+1)=plot(YFC_std.return_period(1:RP100_pos),YFC_std.damage(1:RP100_pos),'or');
    legend_str{end+1} = 'standard storm module';
    legend(legend_handle,legend_str,'Location','SouthEast');
    if sok,saveas(gcf,[fig_dir filesep 'combined_aggregate_DFC_emdat_std.' fig_ext],fig_ext);end
end


%% WISC Synthetic Event Set
% visit 'https://wisc.climate.copernicus.eu/wisc/#/help/products' where you
% download and unpack the synthetic event set footprints (C3S_WISC_EVENT_SET_PART1_0100.tgz to C3S_WISC_EVENT_SET_PART6_0100.tgz)


wisc_dir='D:\Documents_DATA\WISC_data_20170918\Synthetic Event Set'; % or change it to the folder with the unpacked synthetic event set footprints downloaded from WISC



% create hazard out of WISC synthetic event set. save one member each variable/file for better memory matching
[hazard_synthetic1,nc_tryout1]=wisc_event_set_hazard([wisc_dir filesep 'fp_ga3ups_*001.nc'],1,'WISC_SynEventSet_tryout_member1');
clear hazard_synthetic1
[hazard_synthetic2,nc_tryout2]=wisc_event_set_hazard([wisc_dir filesep 'fp_ga3ups_*002.nc'],1,'WISC_SynEventSet_tryout_member2');
clear hazard_synthetic2
[hazard_synthetic3,nc_tryout3]=wisc_event_set_hazard([wisc_dir filesep 'fp_ga3ups_*003.nc'],1,'WISC_SynEventSet_tryout_member3');
clear hazard_synthetic3
[hazard_synthetic4,nc_tryout4]=wisc_event_set_hazard([wisc_dir filesep 'fp_ga3ups_*004.nc'],1,'WISC_SynEventSet_tryout_member4');
clear hazard_synthetic4
[hazard_synthetic5,nc_tryout5]=wisc_event_set_hazard([wisc_dir filesep 'fp_ga3ups_*005.nc'],1,'WISC_SynEventSet_tryout_member5');
clear hazard_synthetic5

hazard_synthetic1 = climada_hazard_load([climada_global.hazards_dir filesep 'WISC_SynEventSet_tryout_member1' '.mat'],1);
hazard_synthetic2 = climada_hazard_load([climada_global.hazards_dir filesep 'WISC_SynEventSet_tryout_member2' '.mat'],1);


% entity=climada_nightlight_global_entity
% entity = climada_entity_load([climada_global.entities_dir filesep 'GLB_0360as_entity.mat'])


entity_combined=[];EDS_era20c=[];EDS_eraint=[];EDS_std_comb=[];  % init
em_data_ET_tot=[];em_data_ET_tot.year=[];em_data_ET_tot.damage=[];em_data_ET_tot.damage_orig=[]; % init
em_data_ST_tot=[];em_data_ST_tot.year=[];em_data_ST_tot.damage=[];em_data_ST_tot.damage_orig=[]; % init


country_names={'Germany','Switzerland'};

for country_i=1:length(country_names)
    
    country_name=country_names{country_i};
    
    % obtain country name and ISO3 code
    [country_name,country_ISO3]=climada_country_name(country_name);
    country_ISO3_title=[country_name ' (' country_ISO3 ')'];
    ISO3_country_name=[country_ISO3 '_' strrep(country_name,' ','')];
    
    fprintf('\n*** processing %s: ***\n',country_ISO3_title);
    
    % create the (default, 10x10km) asset base
    % ----------------------------------------
    
    entity_file=[climada_global.entities_dir filesep 'WISC_' country_ISO3 '_' strrep(country_name,' ','') '_10x10km.mat'];
    if ~exist(entity_file,'file')
        % generate the asset data for requested country on 10x10km resolution
        % ===================================================================
        % (includes the default damage function for storm Europe, too)
        entity=climada_entity_country(country_ISO3);
        entity=climada_assets_encode(entity,hazard_synthetic);
        save(entity_file,'entity',climada_global.save_file_version);
    else
        % load previously generated assets
        entity=climada_entity_load(entity_file,1);
        %save(entity_file,'entity','-v7'); % for Octave compatibility

    end
%     if plot_entity
%         figure;climada_entity_plot(entity); % plot it
%         title(['Assets for ' country_ISO3_title ' (10x10km)'])
%         if sok,saveas(gcf,[fig_dir filesep ISO3_country_name '_assets.' fig_ext],fig_ext);end
%     end % plot_entity
    
    % add to combined entity
    
    % calculate damages for all events for asset base
    % ===============================================
    % EDS = event damage set
    clear EDS % just in case we call this batch code again
    EDS(1)=climada_EDS_calc(entity,hazard_era20c);
    EDS(1).annotation_name='ERA 20c'; % as the name got long
    EDS(2)=climada_EDS_calc(entity,hazard_eraint); % assume same hazard resolution
    EDS(2).annotation_name='ERA int';
    
    for member_i=1:5
        hazard_synthetic = climada_hazard_load([climada_global.hazards_dir filesep 'WISC_SynEventSet_tryout_member' num2str(member_i) '.mat']);
        EDS(member_i+2)=climada_EDS_calc(entity,hazard_synthetic);
        EDS(member_i+2).annotation_name=['synthetic tryout member' num2str(member_i)]; % as the name got long
    end
    
    % plot the exceedance damage frequency curves (DFC)
    % -------------------------------------------------
    
    figure;[~,~,legend_str,legend_handle]=climada_EDS_DFC(EDS);title(['Assets for ' country_ISO3_title])
    if sok,saveas(gcf,[fig_dir filesep ISO3_country_name '_DFC_historic_synthetic.' fig_ext],fig_ext);end
    

end % country_i


climada_global.save_file_version=climada_global_save_file_version; % reset