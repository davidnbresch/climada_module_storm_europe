function severity=winterstorm_compare_severity(compare_scenarios,plot_linear)
% climada
% NAME:
%   winterstorm_compare_severity
% PURPOSE:
%   Compare the hazard severity of different hazard sets
%
%   See also winterstorm_compare
%   See also climada_DFC_compare in core climada
% CALLING SEQUENCE:
%   severity=winterstorm_compare_severity(compare_scenarios)
% EXAMPLE:
%   severity=winterstorm_compare_severity
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   compare_scenarios: if =1, run the severity calculation for single
%       scenarios, too, and add them to the severity comparison
%       =0: omit (default)
%       See HARD CHOICE, TO BE REVISED in code, too
% OUTPUTS:
%   severity: the severity information for all
%       hazard event sets is stored (inspect yourself)
%       For speedup, the routine saves .mat files with each hazard set's
%       results in {module_data_dir}/results
%   plot_linear: if =1, plot linear (default), otherwise loglog (=0)
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141128, initial
%-

severity=[];

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('compare_scenarios','var'),compare_scenarios=0;end
if ~exist('plot_linear','var'),plot_linear=1;end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
%
hazard_set_folder=[module_data_dir filesep 'hazards'];
hazard_set_files={'WS_ECHAM_CTL','WS_ETHC_CTL','WS_GKSS_CTL','WS_ERA40','WS_Europe'}; % last one the blended one
reference_hazard_set=[fileparts(fileparts(climada_global.root_dir)) ...
    filesep 'climada_LOCAL' filesep 'modules' filesep 'WS_Europe_ref' filesep 'WS_reference.mat'];
%
Database_master_table_file=[module_data_dir filesep 'validation' filesep 'Database_master_table.xls'];
%
save_hazard_flag=1; % default=1 for speedup, see winterstorm_scenario_hazard
n_largest_storms=[]; % to limit the storms to compare with to the n largest (defalt=0, all available)
%
% the lines styles
marker_=['*-r';'o-g';'p-b';'.-m';'s- ';'v: ';'d: ';'^: ';'*: ';'o: ';'p--';'.--';'s--';'v--';'d--'];marker_i=1;

severity_i=1;

% for reference hazard
if exist(reference_hazard_set,'file')
    fprintf('reference ');
    load(reference_hazard_set); % load hazard
    
    severity_save_file=[module_data_dir filesep 'results' filesep 'reference_severity.mat'];
    if exist(severity_save_file,'file')
        load(severity_save_file)
        fprintf('loaded\n');
    else
        hazard_severity=winterstorm_severity(hazard); % calculate scenario severity
        save(severity_save_file,'hazard_severity');
    end
    if isempty(severity),clear severity;end
    severity(severity_i)=hazard_severity; % assign
    severity(severity_i).filename='reference';
    severity_i=severity_i+1;
end

% for all hazard event sets
for hazard_i=1:length(hazard_set_files)
    hazard_set_file=[hazard_set_folder filesep hazard_set_files{hazard_i}];
    hazard_set_short=strrep(hazard_set_files{hazard_i},'WS_','');
    fprintf('%s ',hazard_set_short);
    load(hazard_set_file)
    severity_save_file=[module_data_dir filesep 'results' filesep hazard_set_short '_severity.mat'];
    if exist(severity_save_file,'file')
        load(severity_save_file)
        fprintf('loaded\n');
    else
        hazard_severity=winterstorm_severity(hazard); % calculate scenario severity
        save(severity_save_file,'hazard_severity');
    end
    if isempty(severity),clear severity;end
    severity(severity_i)=hazard_severity; % assign
    severity(severity_i).filename=hazard_set_short;
    severity_i=severity_i+1;
end % hazard_i

for hazard_i=1:length(severity)
    [sorted_index,exceedence_freq]=climada_damage_exceedence(severity(hazard_i).index,severity(hazard_i).frequency);
    nonzero_pos     = find(exceedence_freq);
    sorted_index   = sorted_index(nonzero_pos);
    exceedence_freq = exceedence_freq(nonzero_pos);
    return_period   = 1./exceedence_freq;
    legend_str{hazard_i}=strrep(severity(hazard_i).filename,'_',' ');
    if plot_linear
        plot(return_period,sorted_index,marker_(marker_i,:),'markersize',3);hold on
    else
        loglog(return_period,sorted_index,marker_(marker_i,:),'markersize',3);hold on
    end
    marker_i = marker_i+1; if marker_i>length(marker_), marker_i=1; end
end % hazard_i
legend(legend_str)

if compare_scenarios
    % read the master table with storm names and dates
    Database_master_table=climada_spreadsheet_read('no',Database_master_table_file,'table',1);
    
    n_storms=length(Database_master_table.Storm);
    
    if ~isempty(n_largest_storms),n_storms=min(n_storms,n_largest_storms);end
    
    fprintf('comparing with %i scenarios:\n',n_storms);
    
    for storm_i=1:n_storms
        
        storm_data_filename=[module_data_dir filesep 'validation' filesep Database_master_table.Storm{storm_i} '_biasMean.csv'];
        
        % generate the single event hazard set
        hazard=winterstorm_scenario_hazard(storm_data_filename,0,save_hazard_flag);
                
        if ~isempty(hazard)
            
            [~,scenario_name]=fileparts(storm_data_filename);
            scenario_name=strrep(scenario_name,'_',' ');
            scenario_name=strrep(scenario_name,'biasMean','');
            scenario_name=strrep(scenario_name,' ','');
            fprintf('%s ',scenario_name);
            
            severity_save_file=[module_data_dir filesep 'results' filesep scenario_name '_severity.mat'];
            if exist(severity_save_file,'file')
                load(severity_save_file)
                fprintf('loaded\n');
            else
                hazard_severity=winterstorm_severity(hazard); % calculate scenario severity
                save(severity_save_file,'hazard_severity');
            end
            severity(severity_i)=hazard_severity; % assign
            severity(severity_i).filename=scenario_name;
            severity_i=severity_i+1;
            
            hold on;
            XLim = get(get(gcf,'CurrentAxes'),'XLim');
            %YLim = get(get(gcf,'CurrentAxes'),'YLim');
            if plot_linear
                plot(XLim,[hazard_severity.index hazard_severity.index],'-r');
            else
                loglog(XLim,[hazard_severity.index hazard_severity.index],'-r');
            end
            text(XLim(1),hazard_severity.index,scenario_name);
            
        end % ~isempty(hazard)
    end % storm_i
    
end % compare_scenarios

set(gcf,'Color',[1 1 1]);

return
