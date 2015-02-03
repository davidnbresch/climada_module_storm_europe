function hazard=winterstorm_scenario_hazard(storm_data_filename,plot_gust_field,save_hazard_flag)
% climada
% NAME:
%   winterstorm_scenario_hazard
% PURPOSE:
%   read the single storm gust table and create a single storm hazard event
%
%   see winterstorm_validate and winterstorm_compare
% CALLING SEQUENCE:
%   hazard=winterstorm_scenario_hazard(storm_data_filename,save_hazard_flag)
% EXAMPLE:
%   hazard=winterstorm_scenario_hazard(storm_data_filename,save_hazard_flag)
% INPUTS:
%   storm_data_filename: the filename of a single storm (.csv). The file
%   just contains n lines with grid number, gust. Hence we also read the
%   grid table (see )
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   plot_gust_field: if =1, plot the gust field (default=0)
%   save_hazard_flag: if =1, save as single hazard event set (default=0)
%       on subsequent calls, the hazard is read from this file, rather than
%       re-created (hence delete the .mat file to re-create from .csv).
%       If =2, force re-creation of the .mat file, i.e. re-read from the
%       .csv file and create the scenario hazard set again
% OUTPUTS:
%   hazard: a single hazard event set
%       if save_hazard_flag=1, stored as climada hazard event set (with
%       just one event)
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141128, initial
%-

hazard=[]; % init output

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('storm_data_filename','var'),storm_data_filename=[];end
if ~exist('plot_gust_field','var'),plot_gust_field=0;end
if ~exist('save_hazard_flag','var'),save_hazard_flag=0;end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];


% PARAMETERS
%
% define the fil with the grid definiton (grid location, longitude, latitude)
grid_locations_filename=[module_data_dir filesep 'validation' filesep 'grid_locations.csv'];

% prompt for storm_data_filename if not given
if isempty(storm_data_filename) % local GUI
    storm_data_filename=[module_data_dir filesep 'validation' filesep '*.csv'];
    [filename, pathname] = uigetfile(storm_data_filename, 'Select single storm data:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        storm_data_filename=fullfile(pathname,filename);
    end
end

storm_save_filename=strrep(storm_data_filename,'.csv','.mat');

if exist(storm_save_filename,'file') && save_hazard_flag<2
    load(storm_save_filename);
elseif exist(storm_data_filename,'file')
    
    % read the grid on which all gust fields are stored
    RAW = textread(grid_locations_filename,'','delimiter',',','emptyvalue',NaN,'headerlines',1);
    grid.grid_number=RAW(:,1)';
    grid.lon=RAW(:,2)';
    grid.lat=RAW(:,3)';
    
    % read the storm gust field
    RAW = textread(storm_data_filename,'','delimiter',',','emptyvalue',NaN);
    Storm.grid_number=RAW(:,1)';
    Storm.gust=RAW(:,2)';
    
    if length(Storm.grid_number)==length(grid.grid_number) && sum(Storm.grid_number-grid.grid_number)==0
        Storm.lon=grid.lon;
        Storm.lat=grid.lat;
    else
        fprintf('ERROR: storm needs to be gridded, not implemeted yet --> storm skipped\n');
        return
    end
    
    if isfield(Storm,'Longitude')
        
        % create single event hazard set
        
        hazard.comment=sprintf('WSEU event, generated in %s',mfilename);
        hazard.peril_ID='WSEU';
        hazard.date=datestr(now);
        hazard.lat=Storm.lat;
        hazard.lon=Storm.lon;
        hazard.event_count=1;
        hazard.orig_event_flag=1;
        hazard.frequency=1;
        hazard.event_ID=1;
        hazard.orig_years=1;
        hazard.matrix_density=1;
        hazard.centroid_ID=1:length(hazard.lon);
        hazard.filename=storm_data_filename;
        hazard.intensity=Storm.gust;
        
        hazard.intensity=sparse(hazard.intensity);
        
        if save_hazard_flag
            fprintf('saving %s\n',storm_save_filename);
            save(storm_save_filename,'hazard');
        end % save_hazard_flag
        
    end % isfield(Storm,'Longitude')
    
end % exist(storm_data_filename,'file')

if ~isempty(hazard)
    if plot_gust_field
        [~,fN]=fileparts(storm_data_filename);
        fN=strrep(fN,'_',' ');
        climada_color_plot(full(hazard.intensity),hazard.lon,hazard.lat,fN,fN,'','','',0);
    end % plot_gust_field
end

return
