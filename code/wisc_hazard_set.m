function [hazard,nc]=wisc_hazard_set(wisc_file,check_plot,hazard_filename,add_fraction)
% climada WS Europe WISC Copernicus
% MODULE:
%   storm_europe
% NAME:
%   wisc_hazard_set
% PURPOSE:
%   convert WISC storm footprint(s) into climada hazard event set. Visit
%   https://wisc.climate.copernicus.eu/wisc/#/help/products and download
%   C3S_WISC_FOOTPRINT_NETCDF_0100.tgz. This will be your {wisc_dir}. But
%   you can test the code without downloading the data first, see EXAMPLE.
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
%   speed.
%
%   NOTE: if climada_global.save_file_version is set to ='-v7' to allow hazard sets to
%   be read in Octave, too. hazard.fraction is not saved, but re-processed
%   each time the hazard is loaded (saves almost 25% of .mat file size).
%
%   previous call: none
%   next call: wisc_hazard_set_prob, climada_hazard_plot, climada_hazard_ssi, climada_EDS_calc
% CALLING SEQUENCE:
%   [hazard,nc]=wisc_hazard_set(wisc_file,check_plot,hazard_filename,add_fraction)
% EXAMPLE:
%   hazard=wisc_hazard_set('test') % also: wisc_demo
%   hazard=wisc_hazard_set; % create one hazard sewt with all historic events:
%
%   % create hazard sets for ear20c and eraint separately
%   hazard=wisc_hazard_set('{wisc_dir}fp_era20c_*.nc')
%   hazard=wisc_hazard_set('{wisc_dir}fp_eraint_*.nc')
%
%   % in either case, use the following to calculate damages
%   entity=climada_entity_country('DNK');
%   entity=climada_assets_encode(entity,hazard);
%   EDS=climada_EDS_calc(entity,hazard);
%
% INPUTS:
%   wisc_file: the filename (with path) to the WISC footprint data
%       If wisc_file is a regexp (e.g. ='{dir}fp_*.nc'), the code reads all
%       netCDF files according to regexp and stores them in ine hazard
%       event set.
%       If ='test', a TEST windfield (Anatol, Dec 1999) is used (in this
%       special TEST case, the hazard set is stored at the same location as
%       the TEST file, i.e. within the module, not in the hazards folder)
%       If empty, the code assumes ..climada_data/WISC/C3S_WISC_FOOTPRINT_NETCDF_0100/fp_*.nc
%       > promted for if not valid
% OPTIONAL INPUT PARAMETERS:
%   check_plot: if =1, show check plot(s), =0 not (default, except for TEST
%       mode)
%   hazard_filename: the name of the hazard set file, with or without path
%       (in which case path is the default climada_global.hazards_dir).
%       If empty, this is set to a reasonable name
%   add_fraction: if =1, add hazard.fraction. Default =0, since time consuming
%       and not required for WS.
% OUTPUTS:
%   hazard: a climada hazard even set structure, see e.g. climada_tc_hazard_set
%       for a detailed description of all fields. Key fields are:
%       hazard.lon(c): longitude (decimal) of centroid c
%       hazard.lat(c): latitude (decimal) of centroid c
%       hazard.area_km2(c): area of centroid c in km2
%       hazard.on_land(c): =1, if centroid on land, =0 else
%       hazard.lonlat_size: the dimensions of centroids as 2d field. With
%           this, one can convert back to 2d fields by e.g. lon2d=reshape(hazard.lon,hazard.lonlat_size)
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
% David N. Bresch, david.bresch@gmail.com, 20170319, initial
% David N. Bresch, david.bresch@gmail.com, 20170705, climada_entity_country
% David N. Bresch, david.bresch@gmail.com, 20170721, checked with latest WISC data and allow for regexp
% David N. Bresch, david.bresch@gmail.com, 20170730, save for Octave compatibility
% David N. Bresch, david.bresch@gmail.com, 20170801, NaN set to zero
% David N. Bresch, david.bresch@gmail.com, 20171004, checked again, works fine
% Thomas Roeoesli, thomas.roeoesli@usys.ethz.ch, 20171024, added orig_yearset for later use in climada_EDS2YDS
% David N. Bresch, david.bresch@gmail.com, 20171108, double file fp_era20c_1990012515_701_0 excluded
% David N. Bresch, david.bresch@gmail.com, 20180108, see wisc_hazard_set_prob for probabilistic set
% David N. Bresch, david.bresch@gmail.com, 20180109, wisc_file_default better defined
%-

hazard=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('wisc_file','var'),       wisc_file       ='';end
if ~exist('check_plot','var'),      check_plot      = 0;end
if ~exist('hazard_filename','var'), hazard_filename ='';end
if ~exist('add_fraction','var'),    add_fraction    = 0;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% set threshold [m/s] below which we do not store the windfield (for memory reasons)
wind_threshold=15; % default 15 m/s (not too high, to allow to raise a few m/s later in the process if needed, e.g. climate scenarios)
matrix_density=.4; % the average density of the hazard (=#nozero/#all)
%
% whether we create a yearset
create_yearset = true; % default =true
%
% define the TEST wind footprint
%test_wisc_file=[module_data_dir filesep 'raw_windfields' filesep 'fp_3532_1999120318_historic_CF.nc'];
test_wisc_file=[module_data_dir filesep 'raw_windfields' filesep 'fp_eraint_1999122606_507_0.nc']; % better, since same grid as we will use for real
%
% there is one event presetn both in ere20 and eraint for comparison, hence
% we avoid double-counting
exclude_file='fp_era20c_1990012515_701_0.nc';
%
% define the folder and regexp for standard use (i.e. to select all historic events
wisc_file_default=[climada_global.data_dir filesep 'WISC' ...
        filesep 'C3S_WISC_FOOTPRINT_NETCDF_0100' filesep 'fp_era*.nc'];
if isempty(wisc_file),wisc_file=wisc_file_default;end


if ~exist(fileparts(wisc_file),'dir'),wisc_file='';end

%  prompt for wisc_file if not given
if isempty(wisc_file) % local GUI
    wisc_file=[climada_global.data_dir filesep '*.nc'];
    [filename, pathname] = uigetfile(wisc_file, 'Select WISC footprint file:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        wisc_file=fullfile(pathname,filename);
    end
end

TEST_MODE=0;
if strcmpi(wisc_file,'test')
    [wisc_dir,fN,fE]=fileparts(test_wisc_file);
    wisc_files.name=[fN fE];
    % for TEST mode, store the TEST hazard set within the WISC module
    hazard.filename = [wisc_dir filesep fN '.mat'];
    TEST_MODE=1;
else
    wisc_dir=fileparts(wisc_file);
    wisc_files=dir(wisc_file);
    if length(wisc_files)==1
        [~,fN]=fileparts(wisc_files.name);
        hazard.filename = [climada_global.hazards_dir filesep '_' fN '.mat'];
    else
        hazard.filename = [climada_global.hazards_dir filesep 'WISC_eur_WS_hist.mat'];
    end
end

if ~isempty(hazard_filename)
    [fP,fN,fE]=fileparts(hazard_filename);
    if isempty(fP),fP=climada_global.hazards_dir;end
    if isempty(fE),fE='.mat';end
    hazard.filename=[fP filesep fN fE];
end

n_files=length(wisc_files);

if n_files<1
    fprintf(['> visit <a href="https://wisc.climate.copernicus.eu/wisc/#/help/products">'...
        'https://wisc.climate.copernicus.eu/wisc/#/help/products</a>\n',...
        '  and download C3S_WISC_FOOTPRINT_NETCDF_0100.tgz, unzip to find files as\n',...
        '  %s\n'],wisc_file_default);
    return
end

% pre-loop to determine the number of times
% -----------------------------------------
n_events=0;exclude_file_i=0; % init
for file_i=1:n_files
    wisc_file_1=[wisc_dir filesep wisc_files(file_i).name];
    if file_i==1 % take grid info from first file
        fprintf('pre-processing %s',wisc_file_1);
        nc.info = ncinfo(wisc_file_1);
        if ismember('latitude',{nc.info.Variables.Name})
            nc.lat  = ncread(wisc_file_1,'latitude');
            nc.lon  = ncread(wisc_file_1,'longitude');
        else
            nc.lat  = ncread(wisc_file_1,'lat');
            nc.lon  = ncread(wisc_file_1,'lon');
        end
        % construct the 2D grid, if needed
        if size(nc.lat,2)==1,[nc.lat,nc.lon] = meshgrid(nc.lat,nc.lon);end
        
        % add the area [km2] per centroid (111.12 is length of 1 deg in km at equator)
        fprintf('\n> first file, reading lat/lon and calculating centroid area');
        nc.area_km2=abs(nc.lon(2:end,2:end)-nc.lon(1:end-1,1:end-1))...
            .*cos(nc.lat(1:end-1,1:end-1)./180.*pi)*111.12 .* ...
            abs(nc.lat(2:end,2:end)-nc.lat(1:end-1,1:end-1))*111.12;
        % fill last row/column
        nc.area_km2(end+1,:)=nc.area_km2(end,:);
        nc.area_km2(:,end+1)=nc.area_km2(:,end);
    else
        fprintf('pre-processing %s',wisc_files(file_i).name);
    end
    
    % exclude one file that exists twice (for comparison)
    if strcmpi(wisc_files(file_i).name,exclude_file)
        fprintf(' *** excluded (exists only for comparison)');
        exclude_file_i=file_i;
    else
        nc.time = ncread(wisc_file_1,'time');
        n_times = length(nc.time);
        n_events=n_events+n_times; % sum up
    end
    fprintf('\n');
end % file_i

if exclude_file_i>0
    wisc_files(exclude_file_i)=[]; % exclude this file
    n_files=length(wisc_files);
end

n_centroids             = numel(nc.lon);
hazard.lonlat_size      = size(nc.lon);
hazard.lon              = reshape(nc.lon,1,n_centroids);
hazard.lat              = reshape(nc.lat,1,n_centroids);
hazard.area_km2         = reshape(nc.area_km2,1,n_centroids); % to check: hist(hazard.area_km2),xlabel('km^2')

fprintf('> converting %i files with total %i timesteps: ',n_files,n_events);

% the data reading loop
% ---------------------

% allocate the hazard set
hazard.intensity=spalloc(n_events,n_centroids,ceil(n_events*n_centroids*matrix_density));

hazard.yyyy=zeros(1,n_times);
hazard.mm  =zeros(1,n_times);
hazard.dd  =zeros(1,n_times);
hazard.event_ID=zeros(1,n_times);
hazard.orig_event_flag=zeros(1,n_times);

next_i=1; % init
climada_progress2stdout(-1,[],2) % init with mod_step 2 insted of 10
t0 = clock;

for file_i=1:n_files
    
    if length(wisc_files(file_i).name)>2
        wisc_file_1=[wisc_dir filesep wisc_files(file_i).name];
        nc.time = ncread(wisc_file_1,'time');
        n_times = length(nc.time);
        for time_i=1:n_times
            temp_data=ncread(wisc_file_1,'max_wind_gust',[1 1 time_i],[Inf Inf 1]); % only one time slab
            if file_i==1 && time_i==1,nc.wind = temp_data;end % first one stored for checks
            temp_data(temp_data<wind_threshold)=0; % apply winfspeed threshold
            
            hazard.intensity(next_i,:)= sparse(reshape(double(temp_data),1,n_centroids)); % convert 2d -> 1d and sparse
            hazard.event_ID(next_i)   = next_i*100; % we store original events at ID 100,200 ...
            hazard.orig_event_flag(next_i)=1;
            
            % figure date
            temp_str=strrep(wisc_files(file_i).name,'fp_era20c_','');
            temp_str=strrep(               temp_str,'fp_eraint_','');
            hazard.yyyy(next_i) = str2double(temp_str(1:4));
            hazard.mm(  next_i) = str2double(temp_str(5:6));
            hazard.dd(  next_i) = str2double(temp_str(7:8));
            
            climada_progress2stdout(next_i,n_events,10,'events'); % update
            next_i=next_i+1; % explicit, on the safe side
        end % time_i
    end
end % file i
climada_progress2stdout(0) % terminate
t_elapsed_footprints = etime(clock,t0);
fprintf('done, took %3.2f sec. \n',t_elapsed_footprints);

hazard.intensity(isnan(hazard.intensity))=0; % all NaN to zero

if check_plot && TEST_MODE
    fprintf('contour plot ...');
    figure('Name','WISC footprint','Color',[1 1 1]);
    [c,h] = contour(nc.lon,nc.lat,nc.wind);
    clabel(c,h)
    title(strrep(strrep(wisc_files(file_i).name,'.nc',''),'_','\_'));xlabel('max wind gust (m/s)');
    hold on; climada_plot_world_borders(2,'','',1);
    fprintf(' done \n');
end

% Create a yearset to be used in climada_EDS2YDS
if create_yearset
    
    % the beginner does not need to understand whats happening here ;-)
    % see climada_EDS2YDS
    
    year_i=1; % init
    active_year=hazard.yyyy(1); % first year
    event_index=[];event_count=0; % init
    
    fprintf('> generating yearset for %i events: ',n_events);
    climada_progress2stdout    % init, see terminate below
    t0       = clock;
    
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
        
        climada_progress2stdout(footprint_i,n_events,10,'events'); % update
        
    end % footprint_i
    climada_progress2stdout(0) % terminate
    
    % save last year
    hazard.orig_yearset(year_i).yyyy=active_year;
    hazard.orig_yearset(year_i).event_count=event_count;
    hazard.orig_yearset(year_i).event_index=event_index;
    
    t_elapsed_yearset = etime(clock,t0);
    fprintf('done, took %3.2f sec. \n',t_elapsed_yearset);
    
end % create_yearset

% fill/complete hazard structure
hazard.centroid_ID      = 1:n_centroids;
hazard.peril_ID         = 'WS';
hazard.units            = 'm/s';
hazard.date             = datestr(now);
hazard.reference_year   = climada_global.present_reference_year;
hazard.event_count      = n_events;
hazard.orig_event_count = n_events;
hazard.orig_years       = hazard.yyyy(end)-hazard.yyyy(1)+1;
hazard.frequency        = (hazard.event_ID*0+1)/hazard.orig_years;
hazard.comment          = sprintf('WISC WS hazard event set, footprints from %s',wisc_file);
if ~strcmpi(climada_global.save_file_version,'-v7') && add_fraction % gets too big in Octave
    hazard.fraction     = spones(hazard.intensity); % fraction 100%
end
hazard.orig_event_flag  = logical(hazard.orig_event_flag); % convert to logical for indexing
% later to be commented out, for timing
hazard.t_elapsed_footprints = t_elapsed_footprints;
hazard.t_elapsed_yearset    = t_elapsed_yearset;

fprintf('> saving as %s\n',hazard.filename);
save(hazard.filename,'hazard',climada_global.save_file_version);

if check_plot
    if TEST_MODE
        figure('Name','WISC event footprint in CLIMADA','Color',[1 1 1]);
        climada_hazard_plot(hazard,1);title('single wind footprint')
    else
        figure('Name','WISC maximum intensity','Color',[1 1 1]);
        climada_hazard_plot(hazard,0);
    end
end

end % wisc_hazard_set


% OLD code to generate land mask with coastal buffer (not used any more)
% --------
%
% if add_on_land % add a flag to identify centroids on land
%     
%     if exist(hazard_plus.filename,'file')
%         load(hazard_plus.filename)
%         if sum(abs(hazard_plus.lon-hazard.lon))+sum(abs(hazard_plus.lat-hazard.lat))<1000*eps
%             hazard.on_land      = hazard_plus.on_land;
%             hazard.on_landcoast = hazard_plus.on_landcoast;
%             hazard.on_landcoast_comment=hazard_plus.on_landcoast_comment;
%             fprintf('< hazard.on_land and .on_landcoast added from %s\n',hazard_plus.filename);
%         else
%             fprintf('WARNING: hazard.on_land not added (mismatch with grid definition)\n');
%             fprintf('>> delete %s first, re-run %s\n',hazard_plus.filename,mfilename);
%         end
%     else
%         hazard.BoundingBox(1,1)=min(hazard.lon);
%         hazard.BoundingBox(2,1)=max(hazard.lon);
%         hazard.BoundingBox(1,2)=min(hazard.lat);
%         hazard.BoundingBox(2,2)=max(hazard.lat);
%         hazard_rect_x=[hazard.BoundingBox(1,1) hazard.BoundingBox(1,1) hazard.BoundingBox(2,1) hazard.BoundingBox(2,1) hazard.BoundingBox(1,1)];
%         hazard_rect_y=[hazard.BoundingBox(1,2) hazard.BoundingBox(2,2) hazard.BoundingBox(2,2) hazard.BoundingBox(1,2) hazard.BoundingBox(1,2)];
%         fprintf('> adding hazard.on_land (takes time!): ');
%         climada_progress2stdout    % init, see terminate below
%         t0 = clock;
%         hazard.on_land=hazard.lon*0; % init
%         admin0_shapes=climada_shaperead(admin0_shape_file);
%         
%         hazard.country={};
%         for shape_i=1:length(admin0_shapes)
%             % first check for hazard area within country
%             
%             % re-define BoundingBox, as some countries include still their
%             % DOM/TOMs (ie BoundingBox from iriginal shape file, while X
%             % and Y adjusted for climada). Norway (NOR, shape_i=170) is such a case
%             admin0_shapes(shape_i).BoundingBox(1,1)=min(admin0_shapes(shape_i).X);
%             admin0_shapes(shape_i).BoundingBox(2,1)=max(admin0_shapes(shape_i).X);
%             admin0_shapes(shape_i).BoundingBox(1,2)=min(admin0_shapes(shape_i).Y);
%             admin0_shapes(shape_i).BoundingBox(2,2)=max(admin0_shapes(shape_i).Y);
%             
%             country_rect_x=[admin0_shapes(shape_i).BoundingBox(1,1) admin0_shapes(shape_i).BoundingBox(1,1) admin0_shapes(shape_i).BoundingBox(2,1) admin0_shapes(shape_i).BoundingBox(2,1) admin0_shapes(shape_i).BoundingBox(1,1)];
%             country_rect_y=[admin0_shapes(shape_i).BoundingBox(1,2) admin0_shapes(shape_i).BoundingBox(2,2) admin0_shapes(shape_i).BoundingBox(2,2) admin0_shapes(shape_i).BoundingBox(1,2) admin0_shapes(shape_i).BoundingBox(1,2)];
%             hinc=sum(climada_inpolygon(hazard_rect_x,hazard_rect_y,country_rect_x,country_rect_y)); % hazard in country rect
%             cinh=sum(climada_inpolygon(country_rect_x,country_rect_y,hazard_rect_x,hazard_rect_y)); % country in hazard rect
%             if sum(hinc+cinh)>0 % only check for centroids, if BoundingBoxes intersect
%                 fprintf('\n%s\n',admin0_shapes(shape_i).NAME);
%                 country_hit=climada_inpolygon(hazard.lon,hazard.lat,admin0_shapes(shape_i).X,admin0_shapes(shape_i).Y);
%                 hazard.country{end+1}.pos=find(country_hit); % returns index values, saves space
%                 hazard.on_land(country_hit)=1;
%             end
%             climada_progress2stdout(shape_i,length(admin0_shapes),10,'countries'); % update
%         end % shape_i
%         climada_progress2stdout(0) % terminate
%         t_elapsed = etime(clock,t0);
%         fprintf('done, took %3.2f sec. \n',t_elapsed);
%         
%         hazard.on_land=logical(hazard.on_land);
%         
%         % add a boundary zone (coastal waters)
%         hazard.on_landcoast=hazard.on_land;
%         water_points=hazard.on_landcoast == 0;
%         % check with: plot(hazard.lon(water_points),hazard.lat(water_points),'.b','MarkerSize',0.001);hold on
%         distance_km=climada_distance2coast_km(hazard.lon(water_points),hazard.lat(water_points));
%         coastal_points=distance_km<coastal_buffer_km;
%         wp=hazard.on_landcoast(water_points);wp(coastal_points)=1;
%         hazard.on_landcoast(water_points)=wp;
%         % check with: plot(hazard.lon(hazard.on_landcoast),hazard.lat(hazard.on_landcoast),'.g','MarkerSize',0.001);hold on
%         % check with: plot(hazard.lon(hazard.on_land),hazard.lat(hazard.on_land),'.r','MarkerSize',0.0001)
%         
%         fprintf('%2.2i%% water points, thereof %2.2i%% coastal points\n',ceil(sum(water_points)/n_centroids*100),ceil(sum(coastal_points)/sum(water_points)*100));
%         
%         if exist('valid_area','var') % lonmin, lonmax, latmin, latmax
%             fprintf('> restricting to valid_area [%2.2f % 2.2f %2.2f %2.2f]\n',valid_area);
%             valid_area_pos=(hazard.lon>valid_area(1) & hazard.lon<valid_area(2)) & (hazard.lat>valid_area(3) & hazard.lat<valid_area(4));
%             valid_on_landcoast=hazard.on_landcoast(valid_area_pos);
%             hazard.on_landcoast=hazard.on_landcoast*0;
%             hazard.on_landcoast(valid_area_pos)=valid_on_landcoast;
%             hazard.on_landcoast=logical(hazard.on_landcoast);
%             hazard.on_landcoast_comment=sprintf('generated by %s, %s',mfilename,datestr(now));
%         end
% 
%         hazard_plus.lon          = hazard.lon;
%         hazard_plus.lat          = hazard.lat;
%         hazard_plus.lonlat_size  = hazard.lonlat_size;
%         hazard_plus.on_land      = hazard.on_land;
%         hazard_plus.on_landcoast = hazard.on_landcoast;
%         hazard_plus.on_landcoast_comment = hazard.on_landcoast_comment;
%         hazard_plus.area_km2     = hazard.area_km2;
%         
%         % check with: plot(hazard_plus.lon(hazard_plus.on_landcoast),hazard_plus.lat(hazard_plus.on_landcoast),'.g','MarkerSize',0.001);hold on
%         % check with: plot(hazard_plus.lon(hazard_plus.on_land),     hazard_plus.lat(hazard_plus.on_land),'.r','MarkerSize',0.0001)
%         
%         fprintf('> saving area and on_land as as %s\n',hazard_plus.filename);
%         WISC_data_folder=fileparts(hazard_plus.filename);
%         if ~isdir(WISC_data_folder),[fP,fN]=fileparts(WISC_data_folder);mkdir(fP,fN);end % create it
%         save(hazard_plus.filename,'hazard_plus',climada_global.save_file_version);
%         clear hazard_plus
%         
%     end % exist(hazard_plus.filename,'file')
%     
% end % add_on_land