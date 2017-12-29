function [hazard,nc]=wisc_hazard_set(wisc_file,check_plot,hazard_filename,n_prob_events,add_on_land)
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
%   The code is pretty fast, except for detrmining the on_land points.Since
%   this takes long (easiyl an hour, as inpolygon needs to deal with 1.5
%   mio points), the result is saved in a special file in the results
%   folder, named 'WISC_hazard_plus.mat');
%
%   NOTE: if climada_global.save_file_version is set to ='-v7' to allow hazard sets to
%   be read in Octave, too, hazard.fraction is not saved, but re-processed
%   each time the hazard is loaded (saves almost 25% of .mat file size).
%
%   previous call: none
%   next call: climada_EDS_calc or e.g. climada_hazard_plot(hazard,0) or ssi=climada_hazard_ssi(hazard,1)
%   see also: wisc_demo
% CALLING SEQUENCE:
%   [hazard,nc]=wisc_hazard_set(wisc_file,check_plot,hazard_filename,n_prob_events,add_on_land)
% EXAMPLE:
%   hazard=wisc_hazard_set('test') % also: wisc_demo
%   hazard=wisc_hazard_set('test',1,'',4) % smallest probabilistic set
%   % create one hazard sewt with all historic events
%   wisc_dir=[climada_global.data_dir filesep 'WISC' filesep 'C3S_WISC_FOOTPRINT_NETCDF_0100'];
%   hazard=wisc_hazard_set([wisc_dir '__both']) % figures files automatically
%   hazard=wisc_hazard_set([wisc_dir '__both'],0,'',20) % probabilistic set
%
%   % create hazard set for ear20c and eraint separately
%   hazard=wisc_hazard_set('{wisc_dir}fp_era20c_*.nc')
%   hazard=wisc_hazard_set('{wisc_dir}fp_eraint_*.nc')
%
%   % in either case, use the following to calculate damages
%   entity=climada_entity_country('DNK');
%   EDS=climada_EDS_calc(entity,hazard);
%
% INPUTS:
%   wisc_file: the filename (with path) to the WISC footprint data
%       If wisc_file is a regexp (e.g. ='{dir}fn_*.nc'), the code reads all
%       netCDF files according to regexp and stores them in ine hazard
%       event set.
%       If ='test', a TEST windfield (Anatol, Dec 1999) is used (in this
%       special TEST case, the hazard set is stored at the same location as
%       the TEST file, i.e. within the module, not in the hazards folder)
%       If='{dir}__both' (note the two underscores), all files
%       fp_era20c_*.nc and fp_eraint_*.nc in the folder as given in {dir}
%       are processed into one hazard set (in fact, two separate hazard
%       sets are constructed, then merged).
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   check_plot: if =1, show check plot(s), =0 not (default, except for TEST
%       mode)
%   hazard_filename: the name of the hazard set file, with or without path
%       (in which case path is the default climada_global.hazards_dir).
%       If empty, this is set to a reasonable name
%   n_prob_events: default=0, historic set only, if >0, probabilistic
%       events are generated, see also climada_ws_hist2prob
%   add_on_land: if =1, add hazard.on_land (default), =0 not (speeds up)
%       hazard.on_land is used for example for correct ssi calculation (see
%       climada_hazard_ssi). For wisc_file='test', add_on_land is set =0
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
% Thomas R??sli, thomas.roeoesli@usys.ethz.ch, 20171024, added orig_yearset for later use in climada_EDS2YDS
% David N. Bresch, david.bresch@gmail.com, 20171108, double file excluded
% David N. Bresch, david.bresch@gmail.com, 20171229, option '{dir}__both' and fields hazard.area_km2 and hazard.on_land added
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
if ~exist('n_prob_events','var'),   n_prob_events   = 0;end
if ~exist('add_on_land','var'),     add_on_land     = 1;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% set threshold [m/s] below which we do not store the windfield (for memory reasons)
wind_threshold=15; % default 15 m/s (not too high, to allow to raise a few m/s later in the process if needed, e.g. climate scenarios)
matrix_density=.1; % the average density of the hazard (=#nozero/#all)
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
% the default map border file, just in case somebody would like to use something else
admin0_shape_file=climada_global.map_border_file; % only used if add_on_land=1
% the filename where on_land flag is stored (since calculation very time consuming)
hazard_plus.filename=[climada_global.results_dir filesep 'WISC_hazard_plus.mat'];

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

hist_ext='_hist';if n_prob_events>0,hist_ext='';end

TEST_MODE=0;
if strcmpi(wisc_file,'test')
    [wisc_dir,fN,fE]=fileparts(test_wisc_file);
    wisc_files.name=[fN fE];
    % for TEST mode, store the TEST hazard set within the WISC module
    hazard.filename = [wisc_dir filesep fN '.mat'];
    add_on_land=0; % avoid this time-consuming effort
    TEST_MODE=1;
elseif strcmpi(wisc_file(end-5:end),'__both')
    wisc_dir=strrep(wisc_file,'__both','');
    
    hazard.filename=[climada_global.hazards_dir filesep 'WISC_era20c_eur_WS' hist_ext '.mat'];
    if exist(hazard.filename,'file')
        fprintf('loading %s\n',hazard.filename);
        hazard_era20c=climada_hazard_load(hazard.filename);
    else
        hazard_era20c=wisc_hazard_set([wisc_dir filesep 'fp_era20c_*.nc'],0,hazard.filename,n_prob_events,add_on_land);
    end
    
    hazard.filename=[climada_global.hazards_dir filesep 'WISC_eraint_eur_WS' hist_ext '.mat'];
    if exist(hazard.filename,'file')
        fprintf('loading %s\n',hazard.filename);
        hazard_era20c=climada_hazard_load(hazard.filename);
    else
        hazard_eraint=wisc_hazard_set([wisc_dir filesep 'fp_eraint_*.nc'],0,hazard.filename,n_prob_events,add_on_land);
    end
    
    hazard=climada_hazard_merge(hazard_era20c,hazard_eraint,'events');
    hazard=rmfield(hazard,'lonlat_size');
    hazard.filename = [climada_global.hazards_dir filesep 'WISC_eur_WS' hist_ext '.mat'];
    fprintf('> saving combined historic hazard as %s\n',hazard.filename);
    save(hazard.filename,'hazard',climada_global.save_file_version);
    return
else
    wisc_dir=fileparts(wisc_file);
    wisc_files=dir(wisc_file);
    if length(wisc_files)==1
        [~,fN]=fileparts(wisc_files.name);
        hazard.filename = [climada_global.hazards_dir filesep '_' fN '.mat'];
    else
        hazard.filename = [climada_global.hazards_dir filesep '_WISC_eur_WS' hist_ext '.mat'];
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
        '  and download C3S_WISC_FOOTPRINT_NETCDF_0100.tgz\n']);
    return
end

% pre-loop to determine the number of times
% -----------------------------------------
n_events=0;exclude_file_i=0; % init
for file_i=1:n_files
    wisc_file_1=[wisc_dir filesep wisc_files(file_i).name];
    fprintf('pre-processing %s',wisc_file_1);
    if file_i==1 % take grid info from first file
        nc.info = ncinfo(wisc_file_1);
        nc.lat  = ncread(wisc_file_1,'latitude');
        nc.lon  = ncread(wisc_file_1,'longitude');
        % construct the 2D grid, if needed
        if size(nc.lat,2)==1,[nc.lat,nc.lon] = meshgrid(nc.lat,nc.lon);end
        
        % add the area [km2] per centroid (111.12 is length of 1 deg in km at equator)
        fprintf(' (first file, reading lat/lon and calculating centroid area)');
        nc.area_km2=abs(nc.lon(2:end,2:end)-nc.lon(1:end-1,1:end-1))...
            .*cos(nc.lat(1:end-1,1:end-1)./180.*pi)*111.12 .* ...
            abs(nc.lat(2:end,2:end)-nc.lat(1:end-1,1:end-1))*111.12;
        % fill last row/column
        nc.area_km2(end+1,:)=nc.area_km2(end,:);
        nc.area_km2(:,end+1)=nc.area_km2(:,end);
        %hist(reshape(nc.area_km2,1,numel(nc.area_km2))),xlabel('km^2') % to check
        
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

% allocate the hazard set
n_centroids=numel(nc.lon);
hazard.intensity=spalloc(n_events*(n_prob_events+1),n_centroids,ceil(n_events*(n_prob_events+1)*n_centroids*matrix_density));

fprintf('> converting %i files with total %i timesteps: ',n_files,n_events);

% the data reading loop
% ---------------------

hazard.yyyy=zeros(1,(n_prob_events+1)*n_times);
hazard.mm  =zeros(1,(n_prob_events+1)*n_times);
hazard.dd  =zeros(1,(n_prob_events+1)*n_times);
hazard.event_ID=zeros(1,(n_prob_events+1)*n_times);
hazard.orig_event_flag=zeros(1,(n_prob_events+1)*n_times);

next_i=1; % init
climada_progress2stdout    % init, see terminate below
t0       = clock;
for file_i=1:n_files
    if length(wisc_files(file_i).name)>2
        wisc_file_1=[wisc_dir filesep wisc_files(file_i).name];
        nc.time = ncread(wisc_file_1,'time');
        n_times = length(nc.time);
        for time_i=1:n_times
            temp_data=ncread(wisc_file_1,'max_wind_gust',[1 1 time_i],[Inf Inf 1]); % only one time slab
            if file_i==1 && time_i==1,nc.wind = temp_data;end % first one stored for checks
            temp_data(temp_data<wind_threshold)=0;
            if n_prob_events>0
                [intensity_prob,n_prob_events]=climada_ws_hist2prob(double(temp_data));
                i1=1+(next_i-1)*(n_prob_events+1);
                i2=  (next_i  )*(n_prob_events+1);
                hazard.intensity(i1:i2,:) = intensity_prob;
                hazard.event_ID( i1:i2)   = next_i*100+(1:n_prob_events+1)-1;
                hazard.orig_event_flag(i1)=1;
            else
                i1=next_i;i2=i1;
                hazard.intensity(next_i,:)= sparse(reshape(double(temp_data),1,n_centroids));
                hazard.event_ID(next_i)   = next_i*100;
                hazard.orig_event_flag(next_i)=1;
            end
            % figure date
            temp_str=strrep(wisc_files(file_i).name,'fp_era20c_','');
            temp_str=strrep(               temp_str,'fp_eraint_','');
            hazard.yyyy(i1:i2) = str2double(temp_str(1:4));
            hazard.mm(  i1:i2) = str2double(temp_str(5:6));
            hazard.dd(  i1:i2) = str2double(temp_str(7:8));
            
            climada_progress2stdout(next_i,n_events,10,'events'); % update
            next_i=next_i+1; % explicit, on the safe side
        end % time_i
    end
end % file i
climada_progress2stdout(0) % terminate
t_elapsed = etime(clock,t0);
fprintf('done, took %3.2f sec. \n',t_elapsed);

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
    
    t_elapsed = etime(clock,t0);
    fprintf('done, took %3.2f sec. \n',t_elapsed);
    
end % create_yearset

% fill/complete hazard structure
hazard.lonlat_size      = size(nc.lon);
hazard.lon              = reshape(nc.lon,1,n_centroids);
hazard.lat              = reshape(nc.lat,1,n_centroids);
hazard.area_km2         = reshape(nc.area_km2,1,n_centroids); % to check: hist(hazard.area_km2),xlabel('km^2')
hazard.centroid_ID      = 1:n_centroids;
hazard.peril_ID         = 'WS';
hazard.units            = 'm/s';
hazard.date             = datestr(now);
hazard.reference_year   = climada_global.present_reference_year;
hazard.event_count      = n_events*(1+n_prob_events);
hazard.orig_event_count = n_events;
hazard.orig_years       = hazard.yyyy(end)-hazard.yyyy(1)+1;
hazard.frequency        = (hazard.event_ID*0+1)/(hazard.orig_years*(1+n_prob_events));
hazard.comment = sprintf('WISC WS hazard event set, footprints from %s',wisc_file);
if ~strcmpi(climada_global.save_file_version,'-v7') % gets too big in Octave
    hazard.fraction     = spones(hazard.intensity); % fraction 100%
end

if add_on_land % add a flag to identify centroids on land
    
    if exist(hazard_plus.filename,'file')
        load(hazard_plus.filename)
        if sum(abs(hazard_plus.lon-hazard.lon))+sum(abs(hazard_plus.lat-hazard.lat))<1000*eps
            hazard.on_land=hazard_plus.on_land;
        else
            fprintf('WARNING: hazard.on_land not added (mismatch with grid definition)\n');
            fprintf('>> delete %s first, re-run %s\n',hazard_plus.filename,mfilename);
        end
    else
        hazard.BoundingBox(1,1)=min(hazard.lon);
        hazard.BoundingBox(2,1)=max(hazard.lon);
        hazard.BoundingBox(1,2)=min(hazard.lat);
        hazard.BoundingBox(2,2)=max(hazard.lat);
        hazard_rect_x=[hazard.BoundingBox(1,1) hazard.BoundingBox(1,1) hazard.BoundingBox(2,1) hazard.BoundingBox(2,1) hazard.BoundingBox(1,1)];
        hazard_rect_y=[hazard.BoundingBox(1,2) hazard.BoundingBox(2,2) hazard.BoundingBox(2,2) hazard.BoundingBox(1,2) hazard.BoundingBox(1,2)];
        fprintf('> adding hazard.on_land (takes time!): ');
        climada_progress2stdout    % init, see terminate below
        t0 = clock;
        hazard.on_land=hazard.lon*0; % init
        admin0_shapes=climada_shaperead(admin0_shape_file);
        
        for shape_i=1:length(admin0_shapes)
            % first check for hazard area within country
            country_rect_x=[admin0_shapes(shape_i).BoundingBox(1,1) admin0_shapes(shape_i).BoundingBox(1,1) admin0_shapes(shape_i).BoundingBox(2,1) admin0_shapes(shape_i).BoundingBox(2,1) admin0_shapes(shape_i).BoundingBox(1,1)];
            country_rect_y=[admin0_shapes(shape_i).BoundingBox(1,2) admin0_shapes(shape_i).BoundingBox(2,2) admin0_shapes(shape_i).BoundingBox(2,2) admin0_shapes(shape_i).BoundingBox(1,2) admin0_shapes(shape_i).BoundingBox(1,2)];
            hinc=sum(climada_inpolygon(hazard_rect_x,hazard_rect_y,country_rect_x,country_rect_y));
            cinh=sum(climada_inpolygon(country_rect_x,country_rect_y,hazard_rect_x,hazard_rect_y));
            if sum(hinc+cinh)>0 % only check for centroids, if BoundingBoxes intersect
                fprintf('\n%s\n',admin0_shapes(shape_i).NAME);
                country_hit=climada_inpolygon(hazard.lon,hazard.lat,admin0_shapes(shape_i).X,admin0_shapes(shape_i).Y);
                hazard.on_land(country_hit)=1;
            end
            climada_progress2stdout(shape_i,length(admin0_shapes),10,'countries'); % update
        end % shape_i
        climada_progress2stdout(0) % terminate
        t_elapsed = etime(clock,t0);
        fprintf('done, took %3.2f sec. \n',t_elapsed);
        
        % check with
        %plot(hazard.lon(hazard.on_land==1),hazard.lat(hazard.on_land==1),'.r','MarkerSize',0.0001)
        
        hazard_plus.lon         = hazard.lon;
        hazard_plus.lat         = hazard.lat;
        hazard_plus.lonlat_size = hazard.lonlat_size;
        hazard_plus.on_land     = hazard.on_land;
        hazard_plus.area_km2    = hazard.area_km2;
        fprintf('> saving area and on_land as as %s\n',hazard_plus.filename);
        save(hazard_plus.filename,'hazard_plus',climada_global.save_file_version);
        
    end % exist(hazard_plus.filename,'file')
    
end % add_on_land

fprintf('> saving as %s\n',hazard.filename);
save(hazard.filename,'hazard',climada_global.save_file_version);

if check_plot && TEST_MODE
    figure('Name','CLIMADA footprint','Color',[1 1 1]);
    if n_prob_events>1
        subplot(3,3,2),climada_hazard_plot_nogrid(hazard,3);title('North')
        subplot(3,3,4),climada_hazard_plot_nogrid(hazard,4);title('West')
        subplot(3,3,5),climada_hazard_plot_nogrid(hazard,1);title('original')
        subplot(3,3,6),climada_hazard_plot_nogrid(hazard,5);title('East')
        subplot(3,3,8),climada_hazard_plot_nogrid(hazard,2);title('South')
    else
        climada_hazard_plot_nogrid(hazard,1);title('single wind footprint')
    end
end

end % wisc_hazard_set