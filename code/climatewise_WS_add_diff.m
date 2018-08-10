function hazard=climatewise_WS_add_diff(hazard,wsgsmax_diff_file,wsgsmax_diff_var,new_hazard_filename,time_i)
% climada climatewise
% MODULE:
%   storm europe
% NAME:
%   climatewise_WS_add_diff
% PURPOSE:
%   take an existing CLIMADA (WS) hazard event set and add the climate
%   change signal from a netCDF file. The code assumes lon an lat to exist
%   as variables on the file and reads the data from variable
%   WS_wsgsmax_diff_var.
%
%   Stores a new hazard (WS) event set with the name of the original one
%   plus an extension, such as _CC2035
%
%   Usually called from batch-job climatewise_core
%
% CALLING SEQUENCE:
%   hazard=climatewise_WS_add_diff(hazard,wsgsmax_diff_file,wsgsmax_diff_var,new_hazard_filename,time_i);
% EXAMPLE:
%   wsgsmax_diff_file=[climada_global.data_dir filesep 'ClimateWise' filesep 'wsgsmax_diff.nc'];
%   hazard=climatewise_WS_add_diff('WISC_GBR_eur_WS',wsgsmax_diff_file,'wsgsmax_delta_to_baseline');
% INPUTS:
%   hazard: a climada hazard set or a filename (calls climada_hazard_load)
%   wsgsmax_diff_file: a netCDF file with the difference to be added
%       if there is more than one difference field, provide time_i, too
%   wsgsmax_diff_var: the variable name, default is 'wsgsmax_delta_to_baseline'
% OPTIONAL INPUT PARAMETERS:
%   new_hazard_filename: the name of the new hazard set, if not provided,
%       _CC is added. The code does NOT overwrite an existing hazard set,
%       hence delete first, in case you would like to generate the same
%       climate change hazard set again, otherwise the code adds a number
%       to it, such as myhazard_CC02.mat
%   time_i: the timestep (if more than one), default=4 (fors the last time
%       horizon, usually 2045
% OUTPUTS:
%   stores a new hazard set
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20180809, initial
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% init variables (first time called as batch)
if ~exist('hazard','var'),              hazard = [];end
if ~exist('wsgsmax_diff_file','var'),   wsgsmax_diff_file = 'wsgsmax_diff.nc';end
if ~exist('wsgsmax_diff_var','var'),    wsgsmax_diff_var  = 'wsgsmax_delta_to_baseline';end
if ~exist('new_hazard_filename','var'), new_hazard_filename ='';end
if ~exist('time_i','var'),              time_i = 4;end

hazard=climada_hazard_load(hazard);
if isempty(hazard),return;end

% read netCDF file
if exist(wsgsmax_diff_file,'file')
    fprintf('reading lon, lat and %s from %s\n',wsgsmax_diff_var,wsgsmax_diff_file);
    %nc.info          = ncinfo(wsgsmax_diff_file); % info
    nc.time          = ncread(wsgsmax_diff_file,'time');
    nc.lon           = double(ncread(wsgsmax_diff_file,'lon'));
    nc.lat           = double(ncread(wsgsmax_diff_file,'lat'));
    nc.wsgsmax_delta = ncread(wsgsmax_diff_file,wsgsmax_diff_var);
else
    fprintf('ERROR: netCDF file %s not found, aborted\n',wsgsmax_diff_file);
    return
end

% interpolate to centroids
fprintf('interpolating to hazard centroids ...');
wsgsmax_delta=griddata(nc.lon,nc.lat,nc.wsgsmax_delta(:,:,time_i),hazard.lon,hazard.lat,'nearest');
fprintf(' done\n');

% add to each event footprint
n_centroids = size(hazard.intensity,2); % much (!) faster over centroids
fprintf('adding climate change signal to %i centroids (takes some time)\n',n_centroids);
climada_progress2stdout(-1,[],10) % init, see terminate below
for centroid_i=1:n_centroids
    nz_pos=find(hazard.intensity(:,centroid_i));
    hazard.intensity(nz_pos,centroid_i)=hazard.intensity(nz_pos,centroid_i)+wsgsmax_delta(centroid_i);
    climada_progress2stdout(centroid_i,n_centroids,1000,'centroids'); % update
end % event_i
climada_progress2stdout(0) % terminate

hazard.climate_comment=sprintf('%s from %s at time %i',wsgsmax_diff_var,wsgsmax_diff_file,time_i);

if isempty(new_hazard_filename)
    [fP,fN,fE]=fileparts(hazard.filename);
    new_hazard_filename=[fP filesep fN '_CC' fE];
end
CC_count=2; % just to count the hazard set files
while exist(new_hazard_filename,'file') % append to avoid overwriting
    new_hazard_filename=[fP filesep fN '_CC' sprintf('%2.2i',CC_count) fE];CC_count=CC_count+1;
end
fprintf('storing hazard as %s\n',new_hazard_filename);
save(new_hazard_filename,'hazard',climada_global.save_file_version);

end % climatewise_WS_add_diff