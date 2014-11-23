function Database_master_table=winterstorm_validate(plot_validation_chart)
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
%   Storm,Year,Month,Day,Insured_loss,Affected_countries...
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
%   Database_master_table=winterstorm_validate
% EXAMPLE:
%   Database_master_table=winterstorm_validate
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   plot_validation_chart: if =1, plot the validation bar chart (default)
%       =0, do not plot
%       See also plot_gust_field and entity_check_plot in PARAMETERS
%       section to allow for more graphical output.
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
save_hazard_flag=0; % =0

% read the master table with storm names, dates, loss and affected
% countries:
Database_master_table=climada_spreadsheet_read('no',Database_master_table_file,'table',1);
% contains: Storm,Year,Month,Day,Insured_loss,Affected_countries...

% read the grid on which all hust fields are stored
RAW = textread(grid_locations_filename,'','delimiter',',','emptyvalue',NaN,'headerlines',1);
grid.grid_number=RAW(:,1)';
grid.Longitude=RAW(:,2)';
grid.Latitude=RAW(:,3)';

% read the damagefunctions
damagefunctions=climada_damagefunctions_read(damagefunctions_filename);

n_storms=length(Database_master_table.Storm);

fprintf('processing %i storms\n',n_storms);

EDS=[];next_EDS=1; % init

%for storm_i=1:n_storms
for storm_i=1:2
    
    fprintf('%s (reported damage %2.3f bn USD)\n',Database_master_table.Storm{storm_i},Database_master_table.Insured_loss(storm_i)/1e9);
    Modeled_damage=0;
    
    % read the storm gust field
    storm_data_filename=[module_data_dir filesep 'validation' filesep Database_master_table.Storm{storm_i} '_biasMean.csv'];
    if exist(storm_data_filename,'file')
        RAW = textread(storm_data_filename,'','delimiter',',','emptyvalue',NaN);
        Storm.grid_number=RAW(:,1)';
        Storm.gust=RAW(:,2)';
        
        if length(Storm.grid_number)==length(grid.grid_number) && sum(Storm.grid_number-grid.grid_number)==0
            Storm.Longitude=grid.Longitude;
            Storm.Latitude=grid.Latitude;
        else
            fprintf('ERROR: storm needs to be gridded, not implemeted yet --> storm skipped\n');
        end
        
        if isfield(Storm,'Longitude')
            
            if plot_gust_field
                % plot the gust field
                climada_color_plot(Storm.gust,Storm.Longitude,Storm.Latitude,...
                    Database_master_table.Storm{storm_i},...
                    Database_master_table.Storm{storm_i},'','','',0);
            end
            
            % create temporary hazard set
            % ---------------------------
            
            hazard.comment=sprintf('WSEU event, generated in %s',mfilename);
            hazard.peril_ID='WSEU';
            hazard.date=datestr(now);
            hazard.lat=Storm.Latitude;
            hazard.lon=Storm.Longitude;
            hazard.event_count=1;
            hazard.orig_event_flag=1;
            hazard.frequency=1;
            hazard.event_ID=1;
            hazard.orig_years=1;
            hazard.matrix_density=1;
            hazard.centroid_ID=1:length(hazard.lon);
            hazard.filename=storm_data_filename;
            hazard.intensity=Storm.gust;
            
            if save_hazard_flag
                storm_save_filename=strrep(storm_data_filename,'.csv','.mat');
                fprintf('saving %s\n',storm_save_filename);
                save(storm_save_filename,'hazard');
            end
            
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
                    entity.assets.DamageFunID=entity.assets.DamageFunID*0+1;
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
        end
        
    else
        fprintf('ERROR: gust file %s not found, skipped\n',storm_data_filename);
    end
    
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


%bar_chart_table(:,2)=bar_chart_table(:,2)*10; % blow up
%fprintf('WARNING: plotted modeled damages multiplied by 10!\n');
        
if plot_validation_chart
    % show the bar char
    figure ('Name','Event damage comparison')
    set(gcf,'Color',[1 1 1]);
    bar(bar_chart_table,'grouped'); title('Event damage'); legend('Insured','Modeled');
    set(gca,'XTick',1:n_storms);
    set(gca,'XTickLabel',bar_chart_label); % label with storm names
    
    % XTick=get(gca,'XTick'); % following works only for the actual size of the plot
    % for i=1:length(XTick),if XTick(i)>0 && XTick(i)<=length(bar_chart_label),...
    %             Eff_bar_chart_label{i}=bar_chart_label{XTick(i)};end;end
    % set(gca,'XTickLabel',Eff_bar_chart_label); % label with storm names
end

return
