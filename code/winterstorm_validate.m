function Database_master_table=winterstorm_validate(plot_validation_chart,n_largest_storms)
% climada
% NAME:
%   winterstorm_validate
% PURPOSE:
%   Validate the winterstorm module with given historic storms
%
%   See www.europeanwindstorms.org and
%   www.europeanwindstorms.org/cgi-bin/storms/storms.cgi?sort=loss&opt=
%
%   The code first reads the Database_master_table.xls, an Excel file with
%   Storm,Year,Month,Day,Insured_loss,Affected_countries... and the
%   insured value (total value of insurable assets) in Value (used to scale
%   the GDP-based country asset distribution to total insurable asset
%   value).
%
%   Then, it creates the single event hazard set file for each storm, the
%   entity for each of the countries affected by the respective storm and
%   calculates the damage for each country affected by the storm. Finally,
%   it sums up the damage and compares reported and modeled damage (bar
%   chart).
%
%   Note on the insured vaues:
%   A Google search for GDP Germany reveals 3.635 Billionen USD (2013),
%   whereas GDP_entity delivers 3.4951e+12 (extrapolated from 3280529801324
%   USD in 2010 to 2014). Since we need not the GDP, but the insured asset
%   base here,
%
%   See also: winterstorm_calibrate
% CALLING SEQUENCE:
%   Database_master_table=winterstorm_validate(plot_validation_chart,n_largest_storms)
% EXAMPLE:
%   Database_master_table=winterstorm_validate
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   plot_validation_chart: if =1, plot the validation bar chart
%       =0, do not plot (default)
%       See also plot_gust_field and entity_check_plot in PARAMETERS
%       section to allow for more graphical output.
%   n_largest_storms: if not empty, the number of storms treated
%       if empty, all storms are treated (as listed in the
%       Database_master_table, default)
%       Useful for testing, e.g. n_largest_storms=2
% OUTPUTS:
%   Database_master_table: the master table, now with modeled damages and
%       all EDSs (one for each storm for each country)
%       bar_chart_table: the data behind the bar chart (see code)
%       bar_chart_label: the labels of the bar chart (see code)
%   a figure witht the bar chart comparison (make figure window large)
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141122
%-

Database_master_table=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('plot_validation_chart','var'),plot_validation_chart=0;end
if ~exist('n_largest_storms','var'),n_largest_storms=[];end

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% whether we TEST damage functions (as the default damage function seems to
% be too high)
test_damage_function=0; % =0 (no), =1 (DWD) or =2 (EuroTempest), see code below
%
Database_master_table_file=[module_data_dir filesep 'validation' filesep 'Database_master_table.xls'];
grid_locations_filename=[module_data_dir filesep 'validation' filesep 'grid_locations.csv'];
damagefunctions_filename=[module_data_dir filesep 'entities' filesep 'WS_Europe.xls'];
%
unique_DamageFunID=1234; % a kind of unique ID, unlikely to exist
%
% whether we plot the gust field
plot_gust_field=0; % =1
%
% whether we show a check plot of the entity (assets) if newly generated
entity_check_plot=0; % =0
%
% whether we save the single-event (storm) hazard event sets
save_hazard_flag=1; % default=1 (faster)

% read the master table with storm names, dates, loss and affected
% countries:
Database_master_table=climada_spreadsheet_read('no',Database_master_table_file,'table',1);
Value_table=climada_spreadsheet_read('no',Database_master_table_file,'Value',1);
% contains: Storm,Year,Month,Day,Insured_loss,Affected_countries...
Database_master_table.assets=Value_table;
Database_master_table.assets.assets_multiplier=Database_master_table.assets.Value*0; % init

% read the damagefunctions
damagefunctions=climada_damagefunctions_read(damagefunctions_filename);
pos=find(damagefunctions.DamageFunID==1);
damagefunctions.DamageFunID(pos)=damagefunctions.DamageFunID(pos)*0+unique_DamageFunID;

n_storms=length(Database_master_table.Storm);

if ~isempty(n_largest_storms)
    n_storms=min(n_storms,n_largest_storms);
end

fprintf('processing %i storms\n',n_storms);

EDS=[];next_EDS=1; % init

for storm_i=1:n_storms
    
    fprintf('%s (reported damage %2.3f bn USD)\n',Database_master_table.Storm{storm_i},Database_master_table.Insured_loss(storm_i)/1e9);
    Modeled_damage=0;
    
    % read the storm gust field
    storm_data_filename=[module_data_dir filesep 'validation' filesep Database_master_table.Storm{storm_i} '_biasMean.csv'];
    if exist(storm_data_filename,'file')
        
        % generate the single event hazard set
        % ------------------------------------
        hazard=winterstorm_scenario_hazard(storm_data_filename,plot_gust_field,save_hazard_flag);
        
        if ~isempty(hazard)
            
            % figure which country entities we need
            countries=strsplit(Database_master_table.Affected_countries{storm_i},',');
            
            % process affected countries
            % --------------------------
            
            for country_i=1:length(countries)
                
                country_name=strtrim(countries{country_i}); % get rid of insignificant blanks
                
                entity_filename=[climada_global.data_dir filesep 'entities' filesep country_name '_entity.mat'];
                if exist(entity_filename,'file')
                    fprintf('- processing %s',country_name);
                    entity=climada_entity_load(entity_filename);
                else
                    fprintf('- processing %s (creating entity)',country_name);
                    centroids_file     = [climada_global.data_dir filesep 'system'   filesep country_name '_centroids.mat'];
                    entity_file        = [climada_global.data_dir filesep 'entities' filesep country_name '_entity.mat'];
                    entity_future_file = [climada_global.data_dir filesep 'entities' filesep country_name '_entity_future.mat'];
                    
                    % invoke the GDP_entity moduke to generate centroids and entity
                    [centroids,entity,entity_future] = climada_create_GDP_entity(country_name,[],0,1);
                    save(centroids_file,'centroids');
                    save(entity_file,'entity');
                    entity = entity_future; % replace with entity future
                    save(entity_future_file,'entity');
                end
                if entity_check_plot,climada_plot_entity_assets(entity,[],country_name);end
                
                entity_reencoded_file = [climada_global.data_dir filesep 'entities' filesep country_name '_validation_entity.mat'];
                if exist(entity_reencoded_file,'file')
                    load(entity_reencoded_file); % loads entity encoded to validation assets
                else
                    fprintf('\n  re-encoding hazard to the respective centroids\n');
                    assets = climada_assets_encode(entity.assets,hazard);
                    entity=rmfield(entity,'assets');
                    entity.assets=assets; % assign re-encoded assets
                    save(entity_reencoded_file,'entity');
                end
                
                country_pos=strmatch(country_name,Database_master_table.assets.Country);
                if ~isempty(country_pos)
                    if Database_master_table.assets.Value(country_pos)>0
                        assets_Value=sum(entity.assets.Value);
                        assets_multiplier=Database_master_table.assets.Value(country_pos)/assets_Value;
                        entity.assets.Value=entity.assets.Value*assets_multiplier;
                        %fprintf('\n%s assets multiplier %f\n',Database_master_table.assets.Country{country_pos},assets_multiplier);
                        Database_master_table.assets.assets_multiplier(country_pos)=assets_multiplier;
                    else
                        entity.assets.Value=entity.assets.Value*2; % default
                    end
                end
                
                Database_master_table.assets.Value_corrected(country_pos)=sum(entity.assets.Value); % to keep track
                
                if test_damage_function
                    % switch to better damage function
                    if test_damage_function==1
                        % with damage function as derived based on DWD gust data
                        MDR_exp=0.2539;
                        annotation_ext='DWD';
                    elseif test_damage_function==2
                        % with damage function as derived based on Euro Tempest gust data
                        MDR_exp=0.2493;
                        annotation_ext='EuroTempest';
                    end
                    
                    % calculate the damage function
                    entity.assets.DamageFunID=entity.assets.DamageFunID*0+unique_DamageFunID;
                    Intensity=0:2:100;
                    MDR=0.00000001*exp(MDR_exp*Intensity);
                    MDR=min(MDR,1);
                    
                    % switch damage function
                    entity.damagefunctions.DamageFunID=Intensity*0+unique_DamageFunID;
                    entity.damagefunctions.Intensity=Intensity;
                    entity.damagefunctions.MDD=MDR;
                    entity.damagefunctions.PAA=Intensity*0+1;
                    entity.damagefunctions.MDR=MDR;
                    if isfield(entity.damagefunctions,'peril_ID'),entity.damagefunctions=rmfield(entity.damagefunctions,'peril_ID');end
                    
                else
                    % force WS Europe damagefunctions (no surprise)
                    entity=rmfield(entity,'damagefunctions');
                    entity.damagefunctions=damagefunctions;
                    entity.assets.DamageFunID=entity.assets.DamageFunID*0+unique_DamageFunID;
                    annotation_ext='';
                    %fprintf('\n');climada_damagefunctions_map(entity)% for check
                end
                
                % calculate the event loss for country
                % ------------------------------------
                
                EDS_one=climada_EDS_calc(entity,hazard);
                if isempty(EDS)
                    EDS=EDS_one;
                else
                    EDS(next_EDS)=EDS_one;
                end
                EDS(next_EDS).annotation_name=sprintf('%s %s %s',...
                    Database_master_table.Storm{storm_i},country_name,annotation_ext);
                Modeled_damage=Modeled_damage+EDS(next_EDS).ED;
                next_EDS=next_EDS+1;
                
                fprintf(' damage %2.3f bn USD\n',EDS_one.ED/1e9);
                
            end % country_i
        end % ~isempty(hazard)
        
    else
        fprintf('ERROR: gust file %s not found, skipped\n',storm_data_filename);
    end % exist(storm_data_filename,'file')
    
    Database_master_table.Modeled_damage(storm_i)=Modeled_damage;
    fprintf('  modeled %s damage %2.3f bn USD\n',...
        Database_master_table.Storm{storm_i},...
        Database_master_table.Modeled_damage(storm_i)/1e9);
    
    % prepare the data for the bar chart (comparison insured vs modeled)
    bar_chart_table(storm_i,1)=Database_master_table.Insured_loss(storm_i);
    bar_chart_table(storm_i,2)=Database_master_table.Modeled_damage(storm_i);
    
    bar_chart_label{storm_i}=Database_master_table.Storm{storm_i};
    
end % storm_i

Database_master_table.EDS=EDS; % store all EDSs
Database_master_table.bar_chart_table=bar_chart_table; % pass on
Database_master_table.bar_chart_label=bar_chart_label; % pass on

if plot_validation_chart
    % show the bar char
    figure ('Name','Event damage comparison')
    set(gcf,'Color',[1 1 1]);
    bar(bar_chart_table,'grouped','EdgeColor','none'); title('Event damage'); legend('Insured','Modeled');
    set(gca,'XTick',1:length(bar_chart_label));
    set(gca,'XTickLabel',bar_chart_label); % label with storm names
    
    % XTick=get(gca,'XTick'); % following works only for the actual size of the plot
    % for i=1:length(XTick),if XTick(i)>0 && XTick(i)<=length(bar_chart_label),...
    %             Eff_bar_chart_label{i}=bar_chart_label{XTick(i)};end;end
    % set(gca,'XTickLabel',Eff_bar_chart_label); % label with storm names
end

return
