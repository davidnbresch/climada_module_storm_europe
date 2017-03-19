function [hazard,nc]=wisc_hazard_set(wisc_file,check_plot)
% climada WS Europe
% MODULE:
%   storm_europe
% NAME:
%   wisc_hazard_set
% PURPOSE:
%   convert WISC storm footprint(s) into climada hazard even set
%
%   previous call: none
%   next call: climada_EDS_calc
% CALLING SEQUENCE:
%   hazard=wisc_hazard_set(wisc_file)
% EXAMPLE:
%   hazard=wisc_hazard_set('test')
%   entity=climada_nightlight_entity('DNK');
%   EDS=climada_EDS_calc(entity,hazard);
% INPUTS:
%   wisc_file: the filename (with path) to the WISC footprint data
%       If ='test', a TEST windfield (Anatol, Dec 1999) is used
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   check_plot: if =1, show check plot(s), =0 not (default, except for TEST
%       mode)
% OUTPUTS:
%   hazard: a climada hazard even set structure, see e.g.
%       climada_tc_hazard_set for a detailed description of all fields.
%   nc: the content of the netCDF file (for tests, will be disabled later)
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20170319, initial
%-

hazard=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('wisc_file','var'),wisc_file='';end
if ~exist('check_plot','var'),check_plot=0;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% set threshold [m/s] below which we do not store the windfield (for memory reasons)
wind_threshold=15; % default 15 m/s
matrix_density=.1; % the average density of the hazard (=#nozero/#all)
%
% define the TEST wind footprint
test_wisc_file=[module_data_dir filesep 'raw_windfields' filesep 'fp_3532_1999120318_historic_CF.nc'];

% template to prompt for filename if not given
if isempty(wisc_file) % local GUI
    wisc_file=[climada_global.data_dir filesep '*.nc'];
    [filename, pathname] = uigetfile(wisc_file, 'Select WISC footprint file:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        wisc_file=fullfile(pathname,filename);
    end
end

if strcmpi(wisc_file,'test')
    wisc_file=test_wisc_file;
end

[fP,fN]=fileparts(wisc_file);
    
nc.info = ncinfo(wisc_file);
%nc.wind = ncread(wisc_file,'max_wind_gust');
nc.lat  = ncread(wisc_file,'latitude');
nc.lon  = ncread(wisc_file,'longitude');
nc.time = ncread(wisc_file,'time');

n_times = length(nc.time);
n_centroids=numel(nc.lon);

hazard.intensity=spalloc(n_times,n_centroids,ceil(n_times*n_centroids*matrix_density));

fprintf('> converting %i footprints ...',n_times);
for time_i=1:n_times
    temp_data=ncread(wisc_file,'max_wind_gust',[1 1 time_i],[Inf Inf 1]); % only one time slab
    if time_i==1,nc.wind = temp_data;end % first one stored for checks
    temp_data(temp_data<wind_threshold)=0;
    hazard.intensity(time_i,:)=sparse(reshape(double(temp_data),1,n_centroids));
end % time_i
fprintf(' done \n');

if check_plot
    fprintf('contour plot ...');
    figure('Name','WISC footprint','Color',[1 1 1]);
    [c,h] = contour(nc.lon,nc.lat,nc.wind);
    clabel(c,h)
    title(strrep(fN,'_','\_'))
    %climada_plot_world_borders
    fprintf(' done \n');
end

% convert from 2D to list
n_centroids=numel(nc.lon);
n_events=size(nc.wind,3);

% init hazard structure
hazard.lon              = reshape(nc.lon,1,n_centroids);
hazard.lat              = reshape(nc.lat,1,n_centroids);
hazard.centroid_ID      = 1:n_centroids;
hazard.peril_ID         = 'WS';
hazard.units            = 'm/s';
hazard.date             = datestr(now);
hazard.filename         = [fP filesep fN '.mat'];
hazard.reference_year   = climada_global.present_reference_year;
hazard.event_ID         = 1:n_events;
hazard.event_count      = n_events;
hazard.orig_event_flag  = ones(1,n_events);
hazard.orig_event_count = n_events;
hazard.orig_years       = 1;
%nc.info.Variables(1).Attributes(2).Value
%hazard.frequency  = 1/hazard.orig_years;
hazard.frequency  = hazard.event_ID*0+1; % all one
hazard.comment=sprintf('WISC WS hazard event set, footprints from %s',wisc_file);

fprintf('> saving as %s\n',hazard.filename);
save(hazard.filename,'hazard');

end % wisc_hazard_set