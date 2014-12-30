function winterstorm_TEST(check_plots)
% climada
% NAME:
%   winterstorm_TEST
% PURPOSE:
%   TEST European winter storm module
% CALLING SEQUENCE:
%   winterstorm_TEST(check_plot);
% EXAMPLE:
%   winterstorm_TEST(1);
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   check_plots: if =1, show figures to check hazards etc.
%       If =0, skip figures (default)
%       If country_name is set to 'ALL', be careful to set check_plots=1
% OUTPUTS:
%   to stdout and storing an encoded entity to data/entities
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141029
%-

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('check_plots','var'),check_plots=0;end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
WS_entity_file=[module_data_dir filesep 'entities' filesep 'WS_Europe.xls'];
WS_hazard_file=[module_data_dir filesep 'hazards' filesep 'WS_ECHAM_CTL.mat'];

if ~exist(WS_entity_file,'file'),fprintf('ERROR: entity file not found: %s\n',WS_entity_file);return;end
if ~exist(WS_hazard_file,'file'),fprintf('ERROR: hazard set file not found: %s\n',WS_hazard_file);return;end

entity_save_file=strrep(WS_entity_file,'.xls','.mat'); % figure name of the encoded entity file
if ~exist(entity_save_file,'file')
    entity = climada_entity_read(WS_entity_file,WS_hazard_file);
else
    fprintf('%s already exists, delete to re-test\n',entity_save_file);
    load(entity_save_file)
end

load(WS_hazard_file);
EDS=climada_EDS_calc(entity,hazard);
fprintf('WS Europe damage calculation, EL=%f, based on %s\n',EDS.ED,EDS.annotation_name);

if check_plots
    climada_hazard_plot(hazard,0);
end % check_plots

return
