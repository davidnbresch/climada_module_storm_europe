%function climatewise_test
% climada climatewise
% MODULE:
%   storm europe
% NAME:
%   climatewise_test
% PURPOSE:
%   batch job to TEST Climate Wise data etc.
%
%   Input are the bank's assets, stored in single Excel files with one
%   header row, followed by data (all numeric except for postcode_sector,
%   where AB10 1 is the entry, hence this space not to be confused with
%   space as a separator, as it might appear in the raw text here, but all
%   fine within Excel):
%       postcode_sector	latitude	longitude       number of properties    total property value 2016 GBP
%       AB10 1          57.14687255	-2.106558026	17                      3447974
%       AB10 6          57.13707772	-2.122704986	106                     21499132
%       ...
%
%   Please note that for speedup, the hazard once loaded is kept, hence run
%   clear hazard in case you switch to another hazard set. Read the
%   PARAMETERS section carefully anyway (an expert code, not for beginners)
%   Please familiarize yurself with CLIMADA before running this code.
%
% CALLING SEQUENCE:
%   climatewise_test
% EXAMPLE:
%   climatewise_test
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to stdout and figures and to a folder ClimateWise within the CLIMADA
%   results folder
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20170725
% David N. Bresch, david.bresch@gmail.com, 20170729, EM-DAT comparison added
% David N. Bresch, david.bresch@gmail.com, 20180605, new exposure, all NEW
% David N. Bresch, david.bresch@gmail.com, 20180613, absoute and relative
% David N. Bresch, david.bresch@gmail.com, 20180616, paths set to work machine-independent
% David N. Bresch, david.bresch@gmail.com, 20180710, climate scenarios
% David N. Bresch, david.bresch@gmail.com, 20180717, new Excel format, reads replacement_value_gbp instead of total_property_value_2016_GBP
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
%
climate_frequency_screw=0; % NO climate chance scenario
% climate_frequency_screw=1.; % multiplier to the frequency
% climate_intensity_a=1.0;climate_intensity_b=0.2; % intensity=a*intensity+b
%
% whether we process assets (Value, =0) or properties (=1)
process_number_of_properties=1; % default=0
Intensity_threshold_ms=35; % intensity threshold for affected in m/s
%
entity_plot=0; % whether we plot the assets (not needed each time)
upperxlim= 250; % horizontal scale for return period plots
%
% whether we show absolute values (=0) or percentage of Value (=1)
Percentage_Of_Value_Flag=1; % default=0
upperylim=0.25; % vertical scale for return period plots if Percentage_Of_Value_Flag=1
%
% force asset encoding
force_encode=0; % default=0, since automatically detected
%
% define the exposure datasets, one for each bank, just lon/lat/value and
% number_of_properties, will be inserted into a proper entity structure
exposure_folder = [climada_global.data_dir filesep 'ClimateWise'];
exposure_files={
    'barclays_sector_agg.xlsx'
    'hsbc_sector_agg.xlsx'
    'lloyds_sector_agg.xlsx'
    'nationwide_sector_agg.xlsx'
    'rbs_sector_agg.xlsx'
    'santander_sector_agg.xlsx'
    'yorkshire_sector_agg.xlsx'
    };
% local folder to write the figures
fig_dir = [climada_global.results_dir filesep 'ClimateWise'];
if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
fig_ext ='png';
%
hazard_set_file='WISC_GBR_eur_WS'; % fully probabilistic WISC
%hazard_set_file='WISC_eur_WS_hist'; % historic only
%%hazard_set_file='GBR_UnitedKingdom_eur_WS'; % Schwierz et al. combined'best' hazard, same as WS_Europe.mat
%hazard_set_file=[climada_global.modules_dir filesep 'storm_europe/data/hazards' filesep 'WS_Europe.mat']; % Schwierz et al. combined 'best' hazard
%hazard_set_file=[climada_global.modules_dir filesep 'storm_europe/data/hazards' filesep 'WS_ERA40.mat']; % Schwierz et al.
%hazard_set_file=[climada_global.modules_dir filesep 'storm_europe/data/hazards' filesep 'WS_ECHAM_CTL.mat']; % Schwierz et al.
%hazard_set_file=[climada_global.modules_dir filesep 'storm_europe/data/hazards' filesep 'WS_ECHAM_A2.mat']; % Schwierz et al.
%hazard_set_file=[climada_global.modules_dir filesep 'storm_europe/data/hazards' filesep 'WS_ETHC_CTL.mat']; % Schwierz et al.
%hazard_set_file=[climada_global.modules_dir filesep 'storm_europe/data/hazards' filesep 'WS_ETHC_A2.mat']; % Schwierz et al.
%hazard_set_file=[climada_global.modules_dir filesep 'storm_europe/data/hazards' filesep 'WS_GKSS_CTL.mat']; % Schwierz et al.
%hazard_set_file=[climada_global.modules_dir filesep 'storm_europe/data/hazards' filesep 'WS_GKSS_A2.mat']; % Schwierz et al.

% some prep work
if ~exist('hazard','var'),hazard=[];end % init
if isempty(hazard),hazard=climada_hazard_load(hazard_set_file);end
[~,hazard_name]=fileparts(hazard_set_file);
hazard_name=['_' hazard_name]; % prepend _
if strcmpi(hazard_name,'_WISC_GBR_eur_WS')
    hazard_name=''; % default WISC no name
    if Percentage_Of_Value_Flag,upperylim=0.1;end
end

if ~isfield(hazard,'climate_comment') && climate_frequency_screw>0
    hazard.climate_comment=sprintf('CC freq*%2.2f, int*%2.2f+%2.2f',climate_frequency_screw,climate_intensity_a,climate_intensity_b);
    hazard.frequency=hazard.frequency*climate_frequency_screw;
    hazard.climate_frequency_screw=climate_frequency_screw;
    
    if abs(climate_intensity_a-1)+abs(climate_intensity_b)>0
        pos=find(hazard.intensity);
        hazard.intensity(pos)=climate_intensity_a*hazard.intensity(pos)+climate_intensity_b;
        hazard.climate_intensity_a=climate_intensity_a;
        hazard.climate_intensity_b=climate_intensity_b;
    end
    
    fprintf('hazard modified, clear hazard to start with original hazard again\n');
    fprintf('%s\n',hazard.climate_comment);
elseif isfield(hazard,'climate_comment')
    % hazard has been modified, warn the user, if not the same
    if abs(hazard.climate_frequency_screw-climate_frequency_screw)+...
            abs(hazard.climate_intensity_a-climate_intensity_a)+...
            abs(hazard.climate_intensity_b-climate_intensity_b)>0
        fprintf('WARNING: hazard not in line with climate parameters\n');
        fprintf('run clear hazard first\n');
        return
    end
end

Intensity_threshold_ms_str=sprintf('%2.2i',Intensity_threshold_ms);
mode_str='';if process_number_of_properties,mode_str='_properties';upperylim=100;end
Percentage_Of_Value_Flag_str='';if Percentage_Of_Value_Flag,Percentage_Of_Value_Flag_str='_pct';end

entity_template=climada_entity_read('entity_template.xlsx','NOENCODE');
entity_template.assets=rmfield(entity_template.assets,'hazard');
entity_template.assets=rmfield(entity_template.assets,'centroid_index');
entity_template.assets=rmfield(entity_template.assets,'Category_ID');
entity_template.assets=rmfield(entity_template.assets,'Region_ID');
entity_template.assets=rmfield(entity_template.assets,'Value_unit');

entity_template.assets.reference_year=2016;
entity_template.assets.Value_unit{1}='GBP';
entity_template.assets.currency_unit=1e6; % all in mio

% now we start

n_exposures=length(exposure_files);
EDS=[];clear EDS % to re-init
%n_exposures=2; % TEST

if process_number_of_properties
    fprintf('importing %i exposure data sets (processing number_of_properties, threshold %s m/s, %s):\n',n_exposures,Intensity_threshold_ms_str,hazard_name);
else
    fprintf('importing %i exposure data sets (processing Value, %s):\n',n_exposures,hazard_name);
end

for exposure_i=1:n_exposures
    
    exposure_name=strrep(exposure_files{exposure_i},'_sector_agg.xlsx','');
    entity_savefile=[climada_global.entities_dir filesep '_ClimateWise_' exposure_name '.mat'];
    
    if ~exist(entity_savefile,'file')
        fprintf('\nprocessing %s:\n',exposure_name);
        
        exposure_filename=[exposure_folder filesep exposure_files{exposure_i}];
        exposure_data=climada_xlsread('no',exposure_filename);
        
        entity=entity_template;
        entity.assets.filename = exposure_filename;
        entity.assets.lon      = exposure_data.longitude';
        entity.assets.lat      = exposure_data.latitude';
        if isfield(exposure_data,'total_property_value_2016_GBP')
            entity.assets.Value    = exposure_data.total_property_value_2016_GBP'; % until 20180716
        else
            entity.assets.Value    = exposure_data.replacement_value_gbp'; % 20180717
        end
        
        entity.assets.Value=entity.assets.Value/entity.assets.currency_unit;
        
        % complete fields
        entity.assets.Deductible=entity.assets.Value*0;
        entity.assets.Cover=entity.assets.Value;
        entity.assets.DamageFunID=entity.assets.lon*0+1;
        if isfield(exposure_data,'postcode_sector')
            entity.assets.postcode_sector=exposure_data.postcode_sector';
        end
        entity.assets.number_of_properties=exposure_data.number_of_properties';
        entity.assets = climada_assets_complete(entity.assets);
        
        entity=climada_assets_encode(entity,hazard);
        
        save(entity_savefile,'entity');
    else
        fprintf('\nloading %s:\n',exposure_name);
        load(entity_savefile)
    end
    
    % check whether we need to encode the entity's assets:
    if ~strcmpi(entity.assets.hazard.filename,hazard.filename) || force_encode
        encode_entity,entity=climada_assets_encode(entity,hazard);
    end
    
    if process_number_of_properties
        entity.assets.Value=entity.assets.number_of_properties;
        entity.assets.Cover=entity.assets.Value;
        entity.assets.Value_unit=repmat({'properties'},size(entity.assets.Value));
        entity.assets.currency_unit=1; %
        
        entity.damagefunctions.filename='local';
        entity.damagefunctions.Intensity=0:120;
        entity.damagefunctions.MDD=entity.damagefunctions.Intensity*0;
        entity.damagefunctions.PAA=entity.damagefunctions.Intensity*0;
        entity.damagefunctions.DamageFunID=entity.damagefunctions.Intensity*0+1;
        pos=find(entity.damagefunctions.Intensity>=Intensity_threshold_ms);
        entity.damagefunctions.MDD(pos)=1;
        entity.damagefunctions.PAA(pos)=1;
        entity.damagefunctions.peril_ID=repmat({'WS'},size(entity.damagefunctions.Intensity));
        entity.damagefunctions.Intensity_unit=repmat({'m/s'},size(entity.damagefunctions.Intensity));
        entity.damagefunctions.name=repmat({'threshold'},size(entity.damagefunctions.Intensity));
        if isfield(entity.damagefunctions,'datenum'),entity.damagefunctions=rmfield(entity.damagefunctions,'datenum');end
    end % process_number_of_properties
    
    if entity_plot
        figure;climada_entity_plot(entity);
        title(exposure_name);
        saveas(gcf,[fig_dir filesep exposure_name '_exposure' mode_str],fig_ext);
    end % entity_plot
    
    EDS(exposure_i)=climada_EDS_calc(entity,hazard,exposure_name);
    
end % exposure_i

[EDS(end+1),ok]=climada_EDS_combine(EDS);
EDS(end).annotation_name='combined';EDS(end).Value=0;
for EDS_i=1:length(EDS)-1 % sum up exposure
    EDS(end).Value=EDS(end).Value+EDS(EDS_i).Value;
end % EDS_i

figure;climada_EDS_DFC(EDS,[],Percentage_Of_Value_Flag);
xlim([0 upperxlim]);title('Damage exceedance frequency curve')
if Percentage_Of_Value_Flag
    legend('Location','southeast');
    ylim([0 upperylim])
else
    legend('Location','northeast');
end
if process_number_of_properties
    title_str1=['Damage exceedance frequency curve (gust >' Intensity_threshold_ms_str ' m/s)'];
    title(title_str1)
    if isfield(hazard,'climate_comment'),title({title_str1,strrep(hazard.climate_comment,'_',' ')});end
    ylabel('% of properties affected')
    saveas(gcf,[fig_dir filesep  'ClimateWise_EDS' mode_str Percentage_Of_Value_Flag_str '_' Intensity_threshold_ms_str hazard_name],fig_ext);
    return
else
    saveas(gcf,[fig_dir filesep  'ClimateWise_EDS' mode_str Percentage_Of_Value_Flag_str hazard_name],fig_ext);
end


% comparison with UK market portfolio
% -----------------------------------
entity_market=climada_entity_load('GBR_UnitedKingdom_10x10');
entity_market.assets.Value=entity_market.assets.Value*0.75; % convert USD to GBP
V_sum_market_orig=sum(entity_market.assets.Value);
entity_market=climada_assets_encode(entity_market,hazard); % encode to actual hazard

entity_market.assets.currency_unit=entity.assets.currency_unit; % same unit
entity_market.assets.Value=entity_market.assets.Value/entity_market.assets.currency_unit;
entity_market.assets.Cover=entity_market.assets.Value;
V_sum_market=sum(entity_market.assets.Value);V_sum=sum(EDS(end).Value);
entity_market.assets.Value=entity_market.assets.Value/V_sum_market*V_sum; % scale
entity_market.assets.Cover=entity_market.assets.Value;
EDS_market=climada_EDS_calc(entity_market,hazard,'Market portfolio, scaled');
figure;[~,~,legend_str,legend_handle]=climada_EDS_DFC(EDS,EDS_market,Percentage_Of_Value_Flag);
if Percentage_Of_Value_Flag
    legend('Location','southeast');
    ylim([0 upperylim])
else
    legend('Location','northeast');
end
xlim([0 upperxlim]);title('Damage exceedance frequency curve')
saveas(gcf,[fig_dir filesep  'ClimateWise_EDS_comparison' Percentage_Of_Value_Flag_str hazard_name],fig_ext);

% comparison with EM-DAT
% ----------------------
em_data=emdat_read('','GBR','-WS',1,1);
if Percentage_Of_Value_Flag % in percent
    em_data.damage          = em_data.damage          /V_sum_market_orig*100;
    em_data.damage_orig     = em_data.damage_orig     /V_sum_market_orig*100;
    em_data.DFC.damage      = em_data.DFC.damage      /V_sum_market_orig*100;
    em_data.DFC_orig.damage = em_data.DFC_orig.damage /V_sum_market_orig*100;
    legend_location='southeast';
    ylim([0 upperylim])
else
    % scale to market share of combined portfolios
    em_data.damage          = em_data.damage          /V_sum_market*V_sum/entity_market.assets.currency_unit;
    em_data.damage_orig     = em_data.damage_orig     /V_sum_market*V_sum/entity_market.assets.currency_unit;
    em_data.DFC.damage      = em_data.DFC.damage      /V_sum_market*V_sum/entity_market.assets.currency_unit;
    em_data.DFC_orig.damage = em_data.DFC_orig.damage /V_sum_market*V_sum/entity_market.assets.currency_unit;
    legend_location='northeast';
end
[legend_str,legend_handle]=emdat_barplot(em_data,'','','EM-DAT indexed',legend_str,legend_handle,legend_location);
title(sprintf('Damage exceedance curve, total asset value GBP %2.f bn',sum(entity_market.assets.Value)/1e3))
saveas(gcf,[fig_dir filesep  'ClimateWise_EDS_EMDAT_comparison' Percentage_Of_Value_Flag_str hazard_name],fig_ext);



% ------------------------------------------------------------------------
% old code, 20170729:
%
% % global scaling factor to bring damage in line with EM-DAT
% % (CRUDE way of a rough calibration)
% damage_scaling_factor=0.4;
%
% % read Barclay data into climada
% % ------------------------------
% entity=climada_entity_read('ClimateWise_Barclays TEST02.xlsx','GBR_UnitedKingdom_eur_WS');
% entity.assets.currency_unit=1e6; % all in mio
% entity.assets.Value=entity.assets.Value/entity.assets.currency_unit;
% entity.assets.Cover=entity.assets.Value;
% climada_entity_plot(entity) % first plot
% figure;hist(log10(entity.assets.Value));
% title('log10(Value)');xlabel('log10(GBP)');ylabel('# data points (rows in Excel)');set(gcf,'Color',[1 1 1]); % second plot
%
% % damage calculation for Barclays
% % -------------------------------
% EDS=climada_EDS_calc(entity,[],'Barclays'); % damage calculation for Barclays
% EDS.damage=EDS.damage*damage_scaling_factor;EDS.ED=EDS.ED*damage_scaling_factor; % crude
%
% % comparison with UK market portfolio
% % -----------------------------------
% % (scaled to match Barclays 'market share')
% entity_market=climada_entity_load('GBR');
% entity_market.assets.currency_unit=1e6; % all in mio
% entity_market.assets.Value=entity_market.assets.Value/entity_market.assets.currency_unit;
% entity_market.assets.Cover=entity_market.assets.Value;
% V_sum_market=sum(entity_market.assets.Value);V_sum=sum(entity.assets.Value);
% entity_market.assets.Value=entity_market.assets.Value/V_sum_market*V_sum; % scale
% entity_market.assets.Cover=entity_market.assets.Value;
% EDS_market=climada_EDS_calc(entity_market,'GBR_UnitedKingdom_eur_WS','Market portfolio, scaled');
% EDS_market.damage=EDS_market.damage*damage_scaling_factor;EDS_market.ED=EDS_market.ED*damage_scaling_factor; % crude
% figure;[~,~,legend_str,legend_handle]=climada_EDS_DFC(EDS,EDS_market); title('Damage exceedance frequency curve') % third plot
% fprintf('annual expected damage %f mio GBP (market scaled %f)\n',EDS.ED/1e6,EDS_market.ED/1e6)
% fprintf('max damage %f bn GBP (market scaled %f)\n',max(EDS.damage)/1e9,max(EDS_market.damage)/1e9)

% % comparison with EM-DAT
% % ----------------------
% em_data=emdat_read('','GBR','-WS',1,1);
% em_data.damage=em_data.damage/V_sum_market*V_sum/entity.assets.currency_unit;
% em_data.damage_orig=em_data.damage_orig/V_sum_market*V_sum/entity.assets.currency_unit;
% em_data.DFC.damage=em_data.DFC.damage/V_sum_market*V_sum/entity.assets.currency_unit;
% em_data.DFC_orig.damage=em_data.DFC_orig.damage/V_sum_market*V_sum/entity.assets.currency_unit;
% [legend_str,legend_handle]=emdat_barplot(em_data,'','','EM-DAT indexed',legend_str,legend_handle);
% title(sprintf('Damage exceedance curve, total asset value GBP %2.f bn',sum(entity.assets.Value)/1e3))
