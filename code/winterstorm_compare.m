function EDS=winterstorm_compare(entity,compare_damage_functions,compare_hazard_sets,compare_scenarios)
% climada
% NAME:
%   winterstorm_compare
% PURPOSE:
%   Run a winterstorm Europe analysis for a given entity with different
%   damage functions and different hazard event sets.
%
%   compare_damage_functions: those derived based upon DWD and EuroTempest gust
%   data. The code tries to use the matching hazard event set (as as
%   created by country_risk_calc). If not found, the code prompts for a WS
%   hazard event set and re-encodes assets.
%
%   compare_hazard_sets, hazard sets from:
%   Schwierz, C., P. K?llner-Heck, E. Zenklusen Mutter, D. N. Bresch,
%   P.-L.Vidale, M. Wild, C., and Sch?r, 2010: Modelling European winter
%   wind storm losses in current and future climate. Climatic Change (2010)
%   101:485?514, doi: 10.1007/s10584-009-9712-1.
%
%   compare_scenarios: show the modeled scenario loss for the events as
%   defined in the data/validation folder in Database_master_table.xls
%   WARNING: does not make much sense, since the scenario hazard intensity
%   is not calibrated with the hazard intensities in the hazard event sets
%   based on Schwierz et al. Code just kept for reference.
%
%   See also winterstorm_validate and winterstorm_compare_severity
%   See also climada_DFC_compare in core climada
% CALLING SEQUENCE:
%   EDS=winterstorm_compare(entity,compare_damage_functions,compare_hazard_sets,compare_scenarios)
% EXAMPLE:
%   EDS=winterstorm_compare('',1,1) % run damage function and hazard set
%       comparisons
% INPUTS:
%   entity: an encoded entity, see climada_entity_read
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   compare_damage_functions: if =1, compare damage functions (default)
%       run 4 times: damage function as in entity, 2 versions with explicit
%       function and using the standard WS damage function, see
%       damagefunctions_filename in PARAMETERS in code
%       =0: omit this
%   compare_hazard_sets: if =1, compare hazard sets using the standard
%       damage function, see damagefunctions_filename in PARAMETERS in code
%       =0: omit this (default)
%       =2, use the damage function as in entity, not the standard one
%       Note that for each damage calculation, the entity is re-encoded (to
%       be on the safe side)
%   compare_scenarios: if =1, show modeled scenario losses, using standard
%       damage function, see damagefunctions_filename in PARAMETERS in code.
%       See PARAMETERS to define the location and name of the
%       Database_master_table. WARNING: results do not make much sense,
%       since scenario hazard intensities are not calibrated.
%       =0: omit (default)
%       =2, use the damage function as in entity, not the standard one
%       Note that for each damage calculation, the entity is re-encoded (to
%       be on the safe side)
%       Note further that compare_scenarios only does not make much sense,
%       but is permitted, as one might want to overlay scenario results to
%       an existing DFC plot.
% OUTPUTS:
%   EDS: the event damage set(s), see climada_EDS_calc
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141121, ICE initial
% David N. Bresch, david.bresch@gmail.com, 20141128, compare_hazard_sets, compare_scenarios
% David N. Bresch, david.bresch@gmail.com, 20141201, WARNING for compare_scenarios added
%-

EDS=[]; % init output
severity_table=[];

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('entity','var'),entity=[];end
if ~exist('compare_damage_functions','var'),compare_damage_functions=1;end
if ~exist('compare_hazard_sets','var'),compare_hazard_sets=0;end
if ~exist('compare_scenarios','var'),compare_scenarios=0;end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
unique_DamageFunID=1234; % a kind of unique ID, unlikely to exist
%
hazard_set_folder=[module_data_dir filesep 'hazards'];
hazard_set_files={'WS_ECHAM_CTL','WS_ETHC_CTL','WS_GKSS_CTL','WS_ERA40','WS_Europe'}; % last one the blended one
reference_hazard_set=[fileparts(fileparts(climada_global.root_dir)) ...
    filesep 'climada_LOCAL' filesep 'modules' filesep 'WS_Europe_ref' filesep 'WS_reference.mat'];
%
Database_master_table_file=[module_data_dir filesep 'validation' filesep 'Database_master_table.xls'];
%
damagefunctions_filename=[module_data_dir filesep 'entities' filesep 'WS_Europe.xls'];
%
save_hazard_flag=1; % default=1 for speedup, see winterstorm_scenario_hazard
n_largest_storms=[]; % to limit the storms to compare with to the n largest (defalt=0, all available)


% prompt for entity if not given
if isempty(entity),entity=climada_entity_load;end
if isempty(entity),return;end

% We use the WS Europe hazard event set as created by country_risk_calc
%
% Note that one would need to re-encode assets to each hazard prior to
% calling the damage calculation (as climada_EDS_calc assumes matching
% order for speedup), unless one knows that all hazard event sets are valid
% on the exact same centroids (the first n elements in the hazard are
% matching the n locations of the assets, while the n+1:end elements in
% hazard are the ones for the buffer around the country). The call to
% centroids_generate_hazard_sets ensures that, and hence no need for
% re-encoding (force_re_encoding=0). But in case you're in doubt, set
% force_re_encoding=1 and check whether you get the same results (if yes,
% very likely no need for force_re_encoding=1, otherwise keep =1).
% Note further that the code checks re-encoding (issuing a WARNING) if it
% looks like assets and hazard do not match and then sets
% force_re_encoding=1
force_re_encoding=0; % default=0

next_EDS=1; % init

if compare_damage_functions
    
    fprintf('compare damage functions:\n');
    
    WS_hazard_event_set_file=[climada_global.data_dir filesep 'hazards' filesep deblank(entity.assets.filename) '_WS_Europe.mat'];
    
    if ~exist(WS_hazard_event_set_file,'file')
        % propt for WS hazard event set
        WS_hazard_event_set_file=[module_data_dir filesep 'hazards' filesep '*.mat'];
        [filename, pathname] = uigetfile(WS_hazard_event_set_file, 'Select a WS hazard event set:');
        if isequal(filename,0) || isequal(pathname,0)
            return; % cancel
        else
            WS_hazard_event_set_file=fullfile(pathname,filename);
            force_re_encoding=1; % to be on the safe side
        end
    end
    
    load(WS_hazard_event_set_file); % loads hazard
    
    % check for entity to be encoded to the hazard event set's centroids
    % in essence, re-encoce if in doubt (as otherwise results might be useless)
    n_centroids=length(entity.assets.Longitude);
    if size(hazard.lon,1)>1,hazard.lon=hazard.lon';hazard.lat=hazard.lat';end % transpose old hazard event sets
    if length(hazard.lon)>=n_centroids
        if sum(abs(entity.assets.Longitude-hazard.lon(1:n_centroids)))+...
                sum(abs(entity.assets.Latitude-hazard.lat(1:n_centroids)))>0
            fprintf('WARNING: assets might not be encoded to hazard centroids\n');
            force_re_encoding=1;
        end
    else
        fprintf('WARNING: assets not encoded to hazard centroids -> re-encoding\n');
        force_re_encoding=1; % sure to re-encode, not even vector lengths match
    end
    
    EDS=climada_EDS_calc(entity,hazard,'entities'' damage fun',force_re_encoding);
    next_EDS=2;
    
    for test_i=1:3
        
        if test_i==1
            % with damage function as derived based on DWD gust data
            MDR_exp=0.2539;
            annotation_name='DWD';
        elseif test_i==2
            % with damage function as derived based on Euro Tempest gust data
            MDR_exp=0.2493;
            annotation_name='EuroTempest';
        elseif test_i==3
            MDR_exp=[]; % to signal we do not need to calculate
            entity.assets.DamageFunID=entity.assets.DamageFunID*0+1;
            % read WS Europe default damagefunction and replace in entity
            [~,entity]=climada_damagefunctions_read(damagefunctions_filename,entity);
            annotation_name='default';
        end
        
        if ~isempty(MDR_exp)
            % calculate the damage function
            entity.assets.DamageFunID=entity.assets.DamageFunID*0+unique_DamageFunID;
            Intensity=0:2:100;
            MDR=0.00000001*exp(MDR_exp*Intensity);
            MDR=min(MDR,1);
            PAA=Intensity*0+1;
            MDD=MDR./PAA;
            
            entity=rmfield(entity,'damagefunctions'); % remove
            % define damage function:
            entity.damagefunctions.DamageFunID=Intensity*0+unique_DamageFunID;
            entity.damagefunctions.Intensity=Intensity;
            entity.damagefunctions.MDD=MDD;
            entity.damagefunctions.PAA=PAA;
            entity.damagefunctions.MDR=MDR;
            entity.damagefunctions.filename=annotation_name;
        end
        
        % calculate EDS
        EDS(next_EDS)=climada_EDS_calc(entity,hazard,annotation_name);
        next_EDS=next_EDS+1;
        
    end % test_i
    
end % compare_damage_functions

if compare_hazard_sets
    
    fprintf('compare %i hazard event sets:\n',length(hazard_set_files)+1);
    
    if compare_hazard_sets==1
        % make sure we run with default damage function
        [~,entity]=climada_damagefunctions_read(damagefunctions_filename,entity);
        fprintf('NOTE: analysis run with default damage function (%s)\n',damagefunctions_filename);
    elseif compare_hazard_sets==2
        fprintf('NOTE: analysis run with entities'' damage function\n');
    end
    
    if exist(reference_hazard_set,'file')
        fprintf('reference\n');
        load(reference_hazard_set); % load hazard
        if isempty(EDS), clear EDS;end
        EDS(next_EDS)=climada_EDS_calc(entity,hazard,'reference',1); % force re-encoding
        next_EDS=next_EDS+1;
    end
    
    for hazard_i=1:length(hazard_set_files)
        hazard_set_file=[hazard_set_folder filesep hazard_set_files{hazard_i} '.mat'];
        hazard_set_short=strrep(hazard_set_files{hazard_i},'WS_','');
        fprintf('%s\n',hazard_set_short);
        if ~exist(hazard_set_file,'file')
            % try, generate the blended WS_Europe hazard set first
            hazard=winterstorm_blend_hazard_event_sets;
        else
            load(hazard_set_file)
        end
        EDS(next_EDS)=climada_EDS_calc(entity,hazard,hazard_set_short,1); % force re-encoding
        next_EDS=next_EDS+1;
    end % hazard_i
    
    
end % compare_hazard_sets

if compare_hazard_sets || compare_damage_functions
    climada_EDS_DFC(EDS,'',1);
end

if compare_scenarios
        
    % read the master table with storm names and dates
    Database_master_table=climada_spreadsheet_read('no',Database_master_table_file,'table',1);
    
    n_storms=length(Database_master_table.Storm);
    
    if ~isempty(n_largest_storms),n_storms=min(n_storms,n_largest_storms);end
    
    fprintf('comparing with %i storms:\n',n_storms);
    fprintf('WARNING: results do not make much sense, since scenario hazard intensities are not calibrated\n');
    
    if compare_scenarios==1
        % make sure we run with default damage function
        [~,entity]=climada_damagefunctions_read(damagefunctions_filename,entity);
        fprintf('NOTE: analysis run with default damage function (%s)\n',damagefunctions_filename);
    elseif compare_scenarios==2
        fprintf('NOTE: analysis run with entities'' damage function\n');
    end
    
    local_re_encoding=1; % see below
    for storm_i=1:n_storms
        
        storm_data_filename=[module_data_dir filesep 'validation' filesep Database_master_table.Storm{storm_i} '_biasMean.csv'];
        
        % generate the single event hazard set
        hazard=winterstorm_scenario_hazard(storm_data_filename,0,save_hazard_flag);
                
        if ~isempty(hazard)
            
            if local_re_encoding % to avoid re-encoding each time
                entity=climada_assets_encode(entity,hazard);
                local_re_encoding=0;
            end
            
            [~,scenario_name]=fileparts(storm_data_filename);
            scenario_name=strrep(scenario_name,'_',' ');
            scenario_name=strrep(scenario_name,'biasMean','');
            fprintf('%s\n',scenario_name);
            
            % calculate scenario loss
            single_EDS=climada_EDS_calc(entity,hazard,scenario_name,0);
            scenario_damage_pct=single_EDS.ED/single_EDS.Value; % since we show DFC with Percentage_Of_Value_Flag
            
            % add to plot
            hold on;
            XLim = get(get(gcf,'CurrentAxes'),'XLim');
            %YLim = get(get(gcf,'CurrentAxes'),'YLim');
            plot(XLim,[scenario_damage_pct scenario_damage_pct],'-r');
            text(XLim(1),scenario_damage_pct,scenario_name);
            
        end % ~isempty(hazard)
    end % storm_i
    
end % compare_scenarios

title(entity.assets.filename);
hold off

return
