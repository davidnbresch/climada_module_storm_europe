%function climatewise_test
% climada climatewise
% MODULE:
%   storm europe
% NAME:
%   climatewise_test
% PURPOSE:
%   batch job to TEST Climate Wise data etc.
%
%   Process the GBR_*.xlsx exposure data in one batch and the commercial ones
%   ??.xlsx in a second one (see exposure_files in PARAMETERS below)
%
%   for GBR_barclays_sector_agg.xlsx etc:
%   process_number_of_properties=1;
%   Percentage_Of_Value_Flag=1;
%
%   for ??.xlsx  etc:
%   process_number_of_properties=0;
%   Percentage_Of_Value_Flag=0;
%
%   Input are the bank's assets, stored in single Excel files with one
%   header row, followed by data (all numeric except for postcode_sector,
%   where AB10 1 is the entry, hence this space not to be confused with
%   space as a separator, as it might appear in the raw text here, but all
%   fine within Excel):
%       postcode_sector	latitude	longitude       number of properties    replacement_value_gbp
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
% David N. Bresch, david.bresch@gmail.com, 20180720, ??.xlsx files, TC and WS, global
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
%
climate_frequency_screw=0; % NO climate chance scenario
%climate_frequency_screw=1.; % multiplier to the frequency
%climate_intensity_a=1.0;climate_intensity_b=0.2; % intensity=a*intensity+b
%
% whether we process assets (Value, =0) or properties (=1)
process_number_of_properties=0; % default=0
Intensity_threshold_ms_WS=35; % intensity threshold for affected in m/s for WS
Intensity_threshold_ms_TC=55; % intensity threshold for affected in m/s for TC
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
%         'GBR_barclays_sector_agg.xlsx' % GBR_ to indicate GBR only (speedup, must be listed first)
%         'GBR_hsbc_sector_agg.xlsx'
%         'GBR_lloyds_sector_agg.xlsx'
%         'GBR_nationwide_sector_agg.xlsx'
%         'GBR_rbs_sector_agg.xlsx'
%         'GBR_santander_sector_agg.xlsx'
%         'GBR_yorkshire_sector_agg.xlsx'
    '02.xlsx' % follow the global exposures
    '03.xlsx'
    '05.xlsx'
    '08.xlsx'
    '12.xlsx'
    '13.xlsx'
    '14.xlsx'
    '15.xlsx'
    '16.xlsx'
    '23.xlsx'
    '24.xlsx'
    '25.xlsx'
    };
%
% admin0 shape file (to figure country exposure)
admin0_shape_file=climada_global.map_border_file; % as we use the admin0 as in next line as default anyway
%
% local folder to write the figures
fig_dir = [climada_global.results_dir filesep 'ClimateWise'];
if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
fig_ext ='png';
%
%hazard_set_file_WS='WISC_GBR_eur_WS'; % fully probabilistic WISC
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
% --------------
if ~exist('hazard','var'),hazard=[];end
if isempty(hazard),hazard.filename='';end

Intensity_threshold_ms_WS_str=sprintf('%2.2i',Intensity_threshold_ms_WS);
Intensity_threshold_ms_TC_str=sprintf('%2.2i',Intensity_threshold_ms_TC);
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
entity_template.assets.hazard.filename='';

admin0_shapes=climada_shaperead(admin0_shape_file);

% now we start

n_exposures=length(exposure_files);
clear EDS % to re-init
clear EDS_TC % to re-init
%n_exposures=2; % TEST

n_perils=1; % default WS (Europe) only, will be set=2 if TC exposure found

if process_number_of_properties
    fprintf('importing %i exposure data sets (processing number_of_properties, threshold %s/%s m/s):\n',n_exposures,Intensity_threshold_ms_WS_str,Intensity_threshold_ms_TC_str);
else
    fprintf('importing %i exposure data sets (processing Value):\n',n_exposures);
end

for exposure_i=1:n_exposures
    
    exposure_name=strrep(exposure_files{exposure_i},'.xlsx','');
    exposure_name=strrep(exposure_name,'_sector_agg','');
    if length(exposure_name)>2
        exposure_ISO3=exposure_name(1:3); % check if exposure holds the single country name
        [country_name,country_ISO3] = climada_country_name(exposure_ISO3);
        if ~isempty(country_ISO3),exposure_name=exposure_name(5:end);end % without ISO3
    else
        country_ISO3='';
    end
    
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
                
        save(entity_savefile,'entity');
    else
        fprintf('\nloading %s:\n',exposure_name);
        load(entity_savefile)
    end
    
    if entity_plot
        figure;climada_entity_plot(entity);
        title(exposure_name);
        saveas(gcf,[fig_dir filesep exposure_name '_exposure' mode_str],fig_ext);
    end % entity_plot
    
    % figure whether we have one hazard, one country or more complex
    hazard_set_file={};assets_hit={}; % init
    if ~isempty(country_ISO3)
        hazard_set_file{1}=['WISC_' country_ISO3 '_eur_WS']; % single country
        assets_hit{1}.pos=1:length(entity.assets.lon);
    else
        
        % figure which countries and perils to deal with
        
        add_TC=0; % start with Europe only
        for shape_i=1:length(admin0_shapes) % loop over all countries
            % plot(admin0_shapes(shape_i).X,admin0_shapes(shape_i).Y,'-r','LineWidth',2);hold on; axis equal % check plot
            country_hit=climada_inpolygon(entity.assets.lon,entity.assets.lat,admin0_shapes(shape_i).X,admin0_shapes(shape_i).Y);
            if sum(country_hit)>0
                fprintf('%i: %s %s %s\n',shape_i,admin0_shapes(shape_i).NAME,admin0_shapes(shape_i).ADM0_A3,admin0_shapes(shape_i).CONTINENT)
                if strcmpi(admin0_shapes(shape_i).CONTINENT,'Europe')
                    hazard_set_file_test=['WISC_' admin0_shapes(shape_i).ADM0_A3 '_eur_WS'];
                    if exist([climada_global.hazards_dir filesep hazard_set_file_test '.mat'],'file')
                        hazard_set_file{end+1}=['WISC_' admin0_shapes(shape_i).ADM0_A3 '_eur_WS'];
                        assets_hit{end+1}.pos=country_hit;
                    end
                else
                    add_TC=1;n_perils=2;
                end
            end
        end % shape_i    
        
        if add_TC,hazard_set_file{end+1}='GLB_0360as_TC';end
        
    end % figure which countries and perils to deal with
    
    n_hazards=length(hazard_set_file);
    
    for hazard_i=1:n_hazards
        
        [~,fN]=fileparts(hazard.filename);
        if ~strcmp(fN,hazard_set_file{hazard_i}) % not the same as already loaded
            
            fprintf('loading %s',hazard_set_file{hazard_i});
            hazard=climada_hazard_load(hazard_set_file{hazard_i});
            fprintf(' done\n');
            
            if strcmpi(hazard.peril_ID,'WS')
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
            end % CC for WS
        else
            fprintf('hazard already loaded (%s)\n',fN);
        end % loading hazard
        
        % check whether we need to encode the entity's assets:
        if ~strcmpi(entity.assets.hazard.filename,hazard.filename) || force_encode
            entity=climada_assets_encode(entity,hazard);
            if n_hazards==1,save(entity_savefile,'entity');end % if only one hazard, we can save encocded
        end
                
        if process_number_of_properties
            entity.assets.Value=entity.assets.number_of_properties;
            entity.assets.Cover=entity.assets.Value;
            entity.assets.Value_unit=repmat({'properties'},size(entity.assets.Value));
            entity.assets.currency_unit=1;
            
            if strcmpi(hazard.peril_ID,'WS')
                Intensity_threshold_ms=Intensity_threshold_ms_WS;
            else
                Intensity_threshold_ms=Intensity_threshold_ms_TC;
            end
            
            % replace damagefunctions
            entity.damagefunctions.filename='local';
            entity.damagefunctions.Intensity=0:120; % to be on the safe side
            entity.damagefunctions.MDD=entity.damagefunctions.Intensity*0;
            entity.damagefunctions.PAA=entity.damagefunctions.Intensity*0;
            entity.damagefunctions.DamageFunID=entity.damagefunctions.Intensity*0+1;
            pos=find(entity.damagefunctions.Intensity>=Intensity_threshold_ms);
            entity.damagefunctions.MDD(pos)=1;
            entity.damagefunctions.PAA(pos)=1;
            entity.damagefunctions.peril_ID=repmat({hazard.peril_ID},size(entity.damagefunctions.Intensity));
            entity.damagefunctions.Intensity_unit=repmat({'m/s'},size(entity.damagefunctions.Intensity));
            entity.damagefunctions.name=repmat({'threshold'},size(entity.damagefunctions.Intensity));
            if isfield(entity.damagefunctions,'datenum'),entity.damagefunctions=rmfield(entity.damagefunctions,'datenum');end
        end % process_number_of_properties
        
        if n_hazards>1 && strcmp(hazard.peril_ID,'WS') % make sure we do not deal with any WS asset twice (close to country border)
            temp_Value=entity.assets.Value; % copy Values
            entity.assets.Value=entity.assets.Value*0; % set all to zero
            entity.assets.Value(assets_hit{hazard_i}.pos)=temp_Value(assets_hit{hazard_i}.pos); % keep only the ones relevant
        end
        
        EDS_local=climada_EDS_calc(entity,hazard,exposure_name);
        
        if strcmpi(hazard.peril_ID,'WS')
            if hazard_i==1
                EDS(exposure_i)=EDS_local;
            else
                EDS(exposure_i)=climada_EDS_combine(EDS(exposure_i),EDS_local); % add event sets
            end
        elseif strcmpi(hazard.peril_ID,'TC')
            EDS_TC(exposure_i)=EDS_local;
            EDS_TC(exposure_i).annotation_name=exposure_name;
        end
        
    end % hazard_i
    EDS(exposure_i).annotation_name=exposure_name;
    EDS(exposure_i).peril_ID='WS';
    
end % exposure_i

for peril_i=1:n_perils
    
    if peril_i==1
        hazard_name='WS';
        Intensity_threshold_ms_str=Intensity_threshold_ms_WS_str;
    else
        hazard_name='TC';
        EDS_WS=EDS; % to save it
        EDS=EDS_TC;
        Intensity_threshold_ms_str=Intensity_threshold_ms_TC_str;
    end
    
    % remove empty EDSs
    EDS_ok_pos=[];
    for EDS_i=1:length(EDS)
        if ~isempty(EDS(EDS_i).damage),EDS_ok_pos=[EDS_ok_pos EDS_i];end
    end % EDS_i
    EDS=EDS(EDS_ok_pos);
    [EDS(end+1),ok]=climada_EDS_combine(EDS);
    %[EDS_combined,ok]=climada_EDS_combine(EDS);
    %EDS(end+1)=EDS_combined(1); % as sometimes combined contains a 2nd empty set
    EDS(end).annotation_name='combined';EDS(end).Value=0;
    for EDS_i=1:length(EDS)-1 % sum up exposure
        EDS(end).Value=EDS(end).Value+EDS(EDS_i).Value;
    end % EDS_i
    
    figure;climada_EDS_DFC(EDS,[],Percentage_Of_Value_Flag);
    xlim([0 upperxlim]);title([hazard_name ' damage exceedance frequency curve'])
    if Percentage_Of_Value_Flag
        legend('Location','southeast');
        ylim([0 upperylim])
    else
        legend('Location','northeast');
    end
    if process_number_of_properties
        title_str1=[hazard_name ' damage exceedance frequency curve (gust >' Intensity_threshold_ms_str ' m/s)'];
        title(title_str1)
        if isfield(hazard,'climate_comment'),title({title_str1,strrep(hazard.climate_comment,'_',' ')});end
        ylabel('% of properties affected')
        saveas(gcf,[fig_dir filesep  'ClimateWise_EDS' mode_str Percentage_Of_Value_Flag_str '_' Intensity_threshold_ms_str '_' hazard_name],fig_ext);
        return
    else
        saveas(gcf,[fig_dir filesep  'ClimateWise_EDS' mode_str Percentage_Of_Value_Flag_str '_' hazard_name],fig_ext);
    end
    
    if peril_i==1 && ~process_number_of_properties
        
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
        xlim([0 upperxlim]);title([hazard_name ' damage exceedance frequency curve'])
        saveas(gcf,[fig_dir filesep  'ClimateWise_EDS_comparison' Percentage_Of_Value_Flag_str '_' hazard_name],fig_ext);
        
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
        saveas(gcf,[fig_dir filesep  'ClimateWise_EDS_EMDAT_comparison' Percentage_Of_Value_Flag_str '_' hazard_name],fig_ext);
    end
    
end % peril_i

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
