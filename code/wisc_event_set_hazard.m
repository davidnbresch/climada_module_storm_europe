function [hazard,nc]=wisc_event_set_hazard(wisc_file,check_plot,hazard_filename)
% climada WS Europe
% MODULE:
%   storm_europe
% NAME:
%   wisc_event_set_hazard
% PURPOSE:
%   convert WISC synthetic event set into climada hazard event set
%
%   In essence, it reads the variables latitude, longitude, time and
%   max_wind_gust from the netCDF file. If latitude and longitude are 1D,
%   it constructs the 2D grid, otherwise it takes the 2D grid as stored
%   already in latitude and longitude. The code assumes all timesteps to be
%   on the same grid for all times, as it reads it from the first file.
%
%   The code loops over the files twice, first to infer all timesteps, as
%   one netCDF file could contain more than one timestep, then a second
%   time over all files and all timesteps therein to extract the wind
%   speed. To keep the size of the hazard variable within a workable range
%   all information of the wind speed footprint over seas is deleted and only
%   the wind speed footprint over land is kept.
%
%   NOTE: if climada_global.save_file_version is set to ='-v7' to allow hazard sets to
%   be read in Octave, too, hazard.fraction is not saved, but re-processed
%   each time the hazard is loaded (saves almost 25% of .mat file size).
%
%   previous call: none
%   next call: climada_EDS_calc or e.g. climada_hazard_plot(hazard,0)
%   see also: wisc_demo
% CALLING SEQUENCE:
%   hazard=wisc_event_set_hazard(wisc_file)
% EXAMPLE:
%   hazard=wisc_event_set_hazard('test') % also: wisc_demo
%   hazard=wisc_event_set_hazard('{dir}fp_ga3ups_*001.nc')
%   hazard=wisc_event_set_hazard('{dir}fp_ga3ups_*005.nc')
%   entity=climada_entity_country('DNK');
%   EDS=climada_EDS_calc(entity,hazard);
% INPUTS:
%   wisc_file: the filename (with path) to the WISC footprint data
%       If wisc_file is a regexp (e.g. ='{dir}fn_*.nc'), the code reads all
%       netCDF files according to regexp and stores them in ine hazard
%       event set.
%       If ='test', a TEST windfield (Anatol, Dec 1999) is used (in this
%       special TEST case, the hazard set is stored at the same location as
%       the TEST file, i.e. within the module, not in the hazards folder)
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   check_plot: if =1, show check plot(s), =0 not (default, except for TEST
%       mode)
%   hazard_filename: the name of the hazard set file, with or without path
%       (in which case path is the default climada_global.hazards_dir).
%       If empty, this is set to t
% OUTPUTS:
%   hazard: a climada hazard even set structure, see e.g. climada_tc_hazard_set
%       for a detailed description of all fields. Key fields are:
%       hazard.lon(c): longitude (decimal) of centroid c
%       hazard.lat(c): latitude (decimal) of centroid c
%       hazard.intensity(e,c): the wind speed for event e at centroid c
%       hazard.units: the physical units of intensity, here m/s
%       hazard.event_ID(e): ID of event e (usually 1..hazard.event_count)
%       hazard.event_count: the number of events
%       hazard.frequency(e): the frequency of event e = 1/hazard.orig_years
%       hazard.yyyy(e): the year of event e
%       hazard.mm(e): the month of event e
%       hazard.dd(e): the day of event e
%       hazard.comment: a free comment, contains the regexp passed to this function
%       NOTE: for performance reasons, lon and lat are stored as vectors,
%        not as 2D grids any more. Therefore, hazard.intensity(e,c) is a 2D
%        matrix, but hazard.intensity(e,:) is a vector of all centroids for
%        event e and hence can be processed fast in vector-based MATLAB.
%       SPECIAL: hazard.lonlat_size allows to convert back to 2D grid, i.e
%        lon2d=reshape(hazard.lon,hazard.lonlat_size)
%        lat2d=reshape(hazard.lat,hazard.lonlat_size)
%        wind=reshape(hazard.intensity(1,c),hazard.lonlat_size)
%   nc: the content of the netCDF file (for check, info of the first file)
% MODIFICATION HISTORY:
% Thomas R??sli, thomas.roeoesli@usys.ethz.ch, 20170921, initial
% 
%-

hazard=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('wisc_file','var'),wisc_file='';end
if ~exist('check_plot','var'),check_plot=0;end
if ~exist('hazard_filename','var'),hazard_filename='';end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% set threshold [m/s] below which we do not store the windfield (for memory reasons)
wind_threshold=15; % default 15 m/s
matrix_density=.2; % the average density of the hazard (=#nozero/#all)
%
% define the TEST wind footprint
test_wisc_file=[module_data_dir filesep 'raw_windfields' filesep 'fp_ga3ups_198603011200_0022_005.nc'];

% template to prompt for filename if not given
if isempty(wisc_file) % local GUI
    wisc_file=[climada_global.data_dir filesep '*.nc'];
    [filename, pathname] = uigetfile(wisc_file, 'Select WISC event set files:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        wisc_file=fullfile(pathname,filename);
    end
end

if strcmpi(wisc_file,'test')
    [wisc_dir,fN,fE]=fileparts(test_wisc_file);
    wisc_files.name=[fN fE];
    % for TEST mode, store the TEST hazard set within the WISC module
    hazard.filename = [wisc_dir filesep fN '.mat'];
else
    wisc_dir=fileparts(wisc_file);
    wisc_files=dir(wisc_file);
    if length(wisc_files)==1
        [~,fN]=fileparts(wisc_files.name);
        hazard.filename = [climada_global.hazards_dir filesep '_' fN '.mat'];
    else
        hazard.filename = [climada_global.hazards_dir filesep '_WISC_eur_WS.mat'];
    end
end

if ~isempty(hazard_filename)
    [fP,fN,fE]=fileparts(hazard_filename);
    if isempty(fP),fP=climada_global.hazards_dir;end
    hazard.filename=[fP filesep fN fE];
end

n_files=length(wisc_files);

% pre-loop to determine the number of times
% -----------------------------------------
n_events=0; % init
for file_i=1:n_files
    wisc_file_1=[wisc_dir filesep wisc_files(file_i).name];
    fprintf('pre-processing %s\n',wisc_file_1);
    if file_i==1 % take grid info from first file        
        nc.info = ncinfo(wisc_file_1);
        nc.lat  = ncread(wisc_file_1,'lat');
        nc.lon  = ncread(wisc_file_1,'lon');
        % construct the 2D grid, if needed
        if size(nc.lat,2)==1,[nc.lat,nc.lon] = meshgrid(nc.lat,nc.lon);end
    end
    nc.time = ncread(wisc_file_1,'time');
    n_times = length(nc.time);
    n_events=n_events+n_times; % sum up
end % file_i

% cut to land area
country_shapes=climada_shaperead(climada_global.map_border_file,1,1); % reads .mat
land_sea_mask = climada_inshape(nc.lat, nc.lon, country_shapes);

% allocate the hazard set
n_centroids=numel(nc.lon);
hazard.intensity=spalloc(n_events,n_centroids,ceil(n_events*n_centroids*matrix_density));

fprintf('> converting %i files with total %i timesteps: ',n_files,n_events);

% the data reading loop
% ---------------------

next_i=1;
climada_progress2stdout    % init, see terminate below
for file_i=1:n_files
    wisc_file_1=[wisc_dir filesep wisc_files(file_i).name];
    nc.time = ncread(wisc_file_1,'time');
    n_times = length(nc.time);
    for time_i=1:n_times
        temp_data=ncread(wisc_file_1,'wind_speed_of_gust',[1 1],[Inf Inf]); % only one time slab
        if file_i==1 && time_i==1,nc.wind = temp_data;end % first one stored for checks
        temp_data(temp_data<wind_threshold)=0;
        temp_data(~land_sea_mask)=0;
        hazard.intensity(next_i,:)=sparse(reshape(double(temp_data),1,n_centroids));
        % figure date
        temp_str=strrep(wisc_files(file_i).name,'fp_ga3ups_','');
        hazard.yyyy(next_i)=str2double(temp_str(1:4));
        hazard.mm(next_i)  =str2double(temp_str(5:6));
        hazard.dd(next_i)  =str2double(temp_str(7:8));
        climada_progress2stdout(next_i,n_events,10,'events'); % update
        next_i=next_i+1; % explicit, on the safe side
    end % time_i
end % file i
climada_progress2stdout(0) % terminate
fprintf('done \n');

hazard.intensity(isnan(hazard.intensity))=0; % all NaN to zero

if check_plot
    fprintf('contour plot ...');
    figure('Name','WISC footprint','Color',[1 1 1]);
    [c,h] = contour(nc.lon,nc.lat,nc.wind);
    clabel(c,h)
    title(strrep(strrep(wisc_files(file_i).name,'.nc',''),'_','\_'));xlabel('max wind gust (m/s)');
    hold on; climada_plot_world_borders(2,'','',1);
    fprintf(' done \n');
end

% Create a yearset to be used in climada_EDS2YDS
create_yearset = true;
hazard.orig_event_flag  = ones(1,n_events);
if create_yearset
    
    % the beginner does not need to understand whats happening here ;-)
    % see climada_EDS2YDS
    t0       = clock;
    
    year_i=1; % init
    active_year=hazard(1).yyyy(1); % first year
    event_index=[];event_count=0; % init
    
    for footprint_i=1:n_events
        
        if hazard.yyyy(footprint_i)==active_year
            if hazard.orig_event_flag(footprint_i)
                % same year, add if original track
                event_count=event_count+1;
                event_index=[event_index footprint_i];
            end
        else
            % new year, save last year
            hazard.orig_yearset(year_i).yyyy=active_year;
            hazard.orig_yearset(year_i).event_count=event_count;
            hazard.orig_yearset(year_i).event_index=event_index;
            year_i=year_i+1;
            % reset for next year
            active_year=hazard.yyyy(footprint_i);
            if hazard.orig_event_flag(footprint_i)
                % same year, add if original track
                event_count=1;
                event_index=footprint_i;
            end
        end
        
        
    end % footprint_i
    
    % save last year
    hazard.orig_yearset(year_i).yyyy=active_year;
    hazard.orig_yearset(year_i).event_count=event_count;
    hazard.orig_yearset(year_i).event_index=event_index;
    
    t_elapsed = etime(clock,t0);
    msgstr    = sprintf('generating yearset took %3.2f sec',t_elapsed);
    
end % create_yearset



% fill/complete hazard structure
hazard.lonlat_size      = size(nc.lon);
hazard.lon              = reshape(nc.lon,1,n_centroids);
hazard.lat              = reshape(nc.lat,1,n_centroids);
hazard.centroid_ID      = 1:n_centroids;
hazard.peril_ID         = 'WS';
hazard.units            = 'm/s';
hazard.date             = datestr(now);
hazard.reference_year   = climada_global.present_reference_year;
hazard.event_ID         = 1:n_events;
hazard.event_count      = n_events;
% hazard.orig_event_flag  = ones(1,n_events);
hazard.orig_event_count = n_events;
hazard.orig_years       = hazard.yyyy(end)-hazard.yyyy(1)+1;
hazard.frequency        = (hazard.event_ID*0+1)/hazard.orig_years;
hazard.comment = sprintf('WISC Synthetic Event Set, footprints from %s',wisc_file);
if ~strcmpi(climada_global.save_file_version,'-v7') % gets too big in Octave
    hazard.fraction     = spones(hazard.intensity); % fraction 100%
end

fprintf('> saving as %s\n',hazard.filename);
save(hazard.filename,'hazard',climada_global.save_file_version);

end % wisc_hazard_set