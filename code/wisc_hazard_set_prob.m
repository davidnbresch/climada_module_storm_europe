function hazard_info=wisc_hazard_set_prob(country_ISO3,hazard,check_plot,add_fraction)
% climada WS Europe WISC Copernicus
% MODULE:
%   storm_europe
% NAME:
%   wisc_hazard_set_prob
% PURPOSE:
%   add probabilistic events to a WISC hazard event set, see also
%   wisc_hazard_set, as the present code requires wisc_hazard_set to be run
%   (and runs it automatically, if not done so).
%
%   procedure:
%   i)  assign centroids to countries
%       create a file WISC_hazard_info (see OUTPUTS) with most info fields
%   ii) generate probabilistic 'daughter' events
%
%   The code is pretty fast, as it uses matfile to store to disk rather
%   than keeping all in memory (and hence swapping a lot).
%
%   NOTE: if climada_global.save_file_version is set to ='-v7' to allow hazard sets to
%   be read in Octave, hazard.fraction is not saved, but re-processed
%   each time the hazard is loaded (saves almost 25% of .mat file size).
%
%   previous call: wisc_hazard_set (automaticall invoked, if not run before)
%   next call: wisc_hazard_stats, climada_hazard_plot
% CALLING SEQUENCE:
%   hazard_info=wisc_hazard_set_prob(country_ISO3,hazard,check_plot,add_fraction)
% EXAMPLE:
%   hazard_info=wisc_hazard_set_prob('test')
%
%   % then e.g. use the following to calculate damages
%   entity=climada_entity_country('DNK');
%   EDS=climada_EDS_calc(entity,hazard_info);
%
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   country_ISO3: the countries we will generate a probabilistic hazard set
%       for, either single country char ='DNK' or a st ={'DNK','NLD'}.
%       If empty, the setting for country.ISO3 in PARAMETERS is used.
%       Helpful especially to generate a probabilistic set for one single
%       country, e.g. country_ISO3='DNK' is the same as hazard='test'.
%   hazard: the historic hazard set, a climada hazard event set or a
%       filename (if path is missing, the default haazard folder is prepended)
%       If empty, the code searches for a WISC_eur_WS_hist hazard set and
%       runs wisc_hazard_set to generate it if not found
%   check_plot: if =1, show check plot(s), =0 not (default, except for TEST
%       mode)
%   add_fraction: if =1 (default), add hazard.fraction. Set =0 if memory
%       troubles and time issus (fractionnot really required for WS, but
%       you would have to make a local copy of climada_EDS_calc...)
% OUTPUTS:
%   hazard_info: a climada hazard even set structure, see e.g.
%       climada_tc_hazard_set for a detailed description of all fields.
%       This data is returned for the FIRST country. Key fields are:
%       hazard.lon(c): longitude (decimal) of centroid c
%       hazard.lat(c): latitude (decimal) of centroid c
%       hazard.area_km2(c): area of centroid c in km2
%       hazard.lonlat_size: the dimensions of centroids as 2D field. With
%           this, one can convert back to 2D fields by e.g. lon2D=reshape(hazard.lon,hazard.lonlat_size)
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
%        lon2D=reshape(hazard.lon,hazard.lonlat_size)
%        lat2D=reshape(hazard.lat,hazard.lonlat_size)
%        footprint=reshape(hazard.intensity(1,:),hazard.lonlat_size)
%   writes a special file 'WISC_hazard_info.mat' with the core information
%       fields of hazard_info (not intensity) to a WISC subfolder of the
%       climada data folder, a struct with fields:  
%       lon(i) and lat(i) for centroid i as in hazard
%       lonlat_size: the original size of the 2D single footprint array
%       area_km2(i): the area for centroid i in km2
%       shape_i(i): the shape number for centroid i
%       shape.ISO3{j}: the ISO3 code for country j
%       shape.shape_i(j): the shape number for country j, to be used as
%           country_pos=find(hazard_info.shape_i==shape.shape_i(j))
%       comment: a free comment with creation date
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20170108, initial, split from wisc_hazard_set
% David N. Bresch, david.bresch@gmail.com, 20190107, re_create_hazard_info added
%-

hazard_info=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('hazard','var'),          hazard       = [];end
if ~exist('check_plot','var'),      check_plot   =  0;end
if ~exist('add_fraction','var'),    add_fraction =  1;end
if ~exist('country_ISO3','var'),    country_ISO3 = '';end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% the number of probabilistic 'daughter' events per event, please check
% with climada_ws_hist2prob for consistency (used here to pre-allocate memory)
n_prob_events=29; % default =29, see climada_ws_hist2prob
%
% define the countries we will generate a probabilistic hazard set for (at least one)
% country names as ISO3 code, if possible to be alphabetically ordered (ABW..ZWE)
country.ISO3={'AUT','BEL','CHE','DEU','DNK','ESP','FRA','GBR','IRL','LUX','NLD','NOR','POL','PRT','SWE'};
%country.ISO3={'DEU','FRA','GBR'}; % name like 'United Kingdom' or ISO3 code like 'GBR'
%
% set threshold [m/s] below which we do not store the windfield (for memory reasons)
wind_threshold=18; % default 18 m/s (not too high, to allow to raise a few m/s later in the process if needed, e.g. climate scenarios)
matrix_density=.4; % the average density of the hazard (=#nozero/#all)
%
% whether we create a yearset
create_yearset = true; % default =true
%
% the default map border file, just in case somebody would like to use something else
admin0_shape_file=climada_global.map_border_file; % only used if add_on_land=1
% the filename where on_land flag is stored (since calculation very time consuming)
hazard_info.filename=[climada_global.data_dir filesep 'WISC' filesep 'WISC_hazard_info.mat'];
%
% define the climada hazard set with historic footprints (automatically
% generated, if not exising, see wisc_hazard_set)
hazard_filename_hist = [climada_global.hazards_dir filesep 'WISC_eur_WS_hist.mat'];
%
% define hazard filename, XXX will be replaced by country ISO3 code for
% each country in country.ISO3
hazard_filename = [climada_global.hazards_dir filesep 'WISC_XXX_eur_WS.mat'];

if strcmp(country_ISO3,'test')
    country_ISO3={'DNK'};
     fprintf('TEST mode for one single country (%s)\n',country_ISO3{1})
     hazard='';
end

if ~isempty(country_ISO3)
    if ischar(country_ISO3)
        country_ISO3_tmp=country_ISO3;clear country_ISO3
        country_ISO3{1}=country_ISO3_tmp;
    end
    country.ISO3=country_ISO3; % overwrite defaults
end

% load the historic set (if not passed in hazard)
if isempty(hazard),hazard=hazard_filename_hist;end % use default historic set
if ischar(hazard) % a filename passed, try to load
    [fP,fN,fE]=fileparts(hazard);
    if isempty(fP),fP=climada_global.hazards_dir;end
    if isempty(fE),fE='.mat';end
    hazard_filename_hist=[fP filesep fN fE];
    if exist(hazard_filename_hist,'file')
        fprintf('< loading %s\n',hazard_filename_hist)
        load(hazard_filename_hist); % load historic hazard
    else
        fprintf('hazard does not exist, calling wisc_hazard_set to generate it\n')
        hazard=wisc_hazard_set;
    end
end

n_events=size(hazard.intensity,1);
n_events_prob=n_events*(n_prob_events+1); % total number of all events, hist and prob
%n_centroids=size(hazard.intensity,2);
n_countries=length(country.ISO3);

admin0_shapes=climada_shaperead(admin0_shape_file); % read country shapes

if exist(hazard_info.filename,'file') % check for consistency with hazard_info from file
    re_create_hazard_info=0; % assume file is ok
    % add some fields from hazard_info
    load(hazard_info.filename) % contains hazard_info
    if isfield(hazard_info,'lon') && isfield(hazard_info,'lat')
        try
            if sum(abs(hazard_info.lon-hazard.lon))+sum(abs(hazard_info.lat-hazard.lat))<1000*eps % check for same grid
                % grid ok, now check for requested countries being available in hazard_info
                shape_is=zeros(1,n_countries);
                for shape_i=1:length(admin0_shapes) % convert ISO3 to shape index, as faster
                    country_contains=contains(country.ISO3,admin0_shapes(shape_i).ADM0_A3);
                    if sum(contains(country.ISO3,admin0_shapes(shape_i).ADM0_A3)) % country within set of requestedones
                        shape_is(country_contains) = shape_i;
                    end
                end % shape_i
                if sum(ismember(shape_is,hazard_info.shape.shape_i))==length(shape_is)
                    fprintf('< hazard_info from %s\n',hazard_info.filename);
                else
                    re_create_hazard_info=1;
                    fprintf('WARNING: hazard_info mismatch with country list, re-created\n');
                end
            end
        catch
            re_create_hazard_info=1;
            fprintf('WARNING: hazard_info corrupted, re-created\n');
        end
    else
        re_create_hazard_info=1;
        fprintf('WARNING: hazard_info mismatch with grid definition, re-created\n');
    end
    if re_create_hazard_info
        hazard_info=rmfield(hazard_info,'shape_i');
        hazard_info=rmfield(hazard_info,'shape');
    end
end

if ~isfield(hazard_info,'shape_i') % (re)create hazard_info
    
    fprintf('assigning centroids for %i countries:\n',n_countries);
    
    % init hazard_info
    hazard_info.lon     = hazard.lon;
    hazard_info.lat     = hazard.lat;
    hazard_info.shape_i = hazard.lon*0; % init
    if isfield(hazard,'lonlat_size'),hazard_info.lonlat_size  = hazard.lonlat_size;end
    if ~isfield(hazard,'area_km2')
        fprintf('Warning: field hazard.area_km2 missing, consider re-running wisc_hazard_set\n');
    else
        hazard_info.area_km2     = hazard.area_km2;
    end
    
    t0=clock;
    for shape_i=1:length(admin0_shapes)
        
        country_contains=contains(country.ISO3,admin0_shapes(shape_i).ADM0_A3);
        if sum(country_contains) % country within set of requestedones
            country_i=find(country_contains);
            
            fprintf('  centroids for %s\n',admin0_shapes(shape_i).ADM0_A3);
            country_hit=climada_inpolygon(hazard_info.lon,hazard_info.lat,admin0_shapes(shape_i).X,admin0_shapes(shape_i).Y);
            hazard_info.shape_i(country_hit)     = shape_i;
            hazard_info.shape.ISO3{country_i}    = admin0_shapes(shape_i).ADM0_A3;
            hazard_info.shape.shape_i(country_i) = shape_i;
            
        end
    end % shape_i
    t_elapsed_shape_i = etime(clock,t0);
    fprintf('done, took %3.2f sec. \n',t_elapsed_shape_i);
    
    fprintf('> saving hazard info to %s\n',hazard_info.filename);
    WISC_data_folder=fileparts(hazard_info.filename);
    if ~isdir(WISC_data_folder),[fP,fN]=fileparts(WISC_data_folder);mkdir(fP,fN);end % create it
    hazard_info.comment=sprintf('created by %s at %s',mfilename,datestr(now)); % to have date as precise as possible
    save(hazard_info.filename,'hazard_info',climada_global.save_file_version);
    
end % exist(hazard_info.filename,'file')

clear admin0_shapes

% figure idex of centroids for each country (for speedup)
if isfield(country,'pos'),country=rmfield(country,'pos');end
for country_i=1:n_countries
    iii=[];
    for si=1:length(hazard_info.shape.shape_i)
        % find country.ISO3 in hazard_info.shape.ISO3
        if strcmpi(hazard_info.shape.ISO3{si},country.ISO3{country_i}),iii=si;end
    end %si
    country.pos{country_i}.pos=find(hazard_info.shape_i==hazard_info.shape.shape_i(iii)); % returns index values, saves space
    % e.g.: hazard.shape_i(country.pos{country_i}.pos)) % all centroids within country i
end % country_i

% init the country hazards as matfiles (for memory reasons)
[fP,fN,fE]=fileparts(hazard_filename);
for country_i=1:n_countries
    fName=strrep(fN,'XXX',country.ISO3{country_i});
    country.hazard_filename{country_i}=[fP filesep fName fE];
    country.matfile{country_i}=[climada_global.results_dir filesep fName fE];
    %hazards{country_i}=matfile(country.matfile{country_i}); % allocate matfile
    n_country_centroids=length(country.pos{country_i}.pos);
    hazards{country_i}.intensity=spalloc(n_events_prob,n_country_centroids,...
        ceil(n_events_prob*n_country_centroids*matrix_density));
end % country_i

% copy from hazard
hazard_info.peril_ID         = hazard.peril_ID;
hazard_info.units            = hazard.units;
hazard_info.reference_year   = hazard.reference_year;
hazard_info.orig_event_count = hazard.orig_event_count;
hazard_info.orig_years       = hazard.orig_years;
hazard_info.centroid_ID      = hazard.centroid_ID;
% allocate
hazard_info.event_ID         = zeros(1,n_events_prob);
hazard_info.yyyy             = zeros(1,n_events_prob);
hazard_info.mm               = zeros(1,n_events_prob);
hazard_info.dd               = zeros(1,n_events_prob);
hazard_info.orig_event_flag  = zeros(1,n_events_prob);
hazard_info.comment          = sprintf('WISC probabilistic WS hazard event set, historic footprints from %s',hazard.filename);
hazard_yyyy=hazard.yyyy;hazard_mm=hazard.mm;hazard_dd=hazard.dd;
hazard_lonlat_size=hazard_info.lonlat_size;
hazard_intensity=hazard.intensity; % better for sparse

% get rid of as many fields in hazard as possible (to save memory)
% hazard=rmfield(hazard,'centroid_ID');
% hazard=rmfield(hazard,'event_ID');
clear hazard;

fprintf('> converting %i events for %i countries: ',n_events,n_countries);

climada_progress2stdout(-1,[],2) % init with mod_step 2 insted of 10

t0 = clock;
for event_i=1:n_events % for fast TEST, use event_i=55:56 for DNK
    
    intensity_hist1D=hazard_intensity(event_i,:); % make full
    intensity2D=reshape(full(intensity_hist1D),hazard_lonlat_size); % covert to 2D
    intensity_prob=climada_ws_hist2prob(intensity2D,'',n_prob_events); % input is 2D, output is 1D
    intensity_prob(intensity_prob<wind_threshold)=0; % set low winds to zero
    
    i1=1+(event_i-1)*(n_prob_events+1);
    i2=  (event_i  )*(n_prob_events+1);
    
    for country_i=1:n_countries
        hazards{country_i}.intensity(i1     ,:) = intensity_hist1D(1,country.pos{country_i}.pos); % original event
        hazards{country_i}.intensity(i1+1:i2,:) = sparse(intensity_prob(:,country.pos{country_i}.pos));
    end % country_i
    
    climada_progress2stdout(event_i,n_events_prob,10,'events'); % update
end % event_i
t_elapsed_footprints = etime(clock,t0);
climada_progress2stdout(0) % terminate
fprintf('done, took %3.2f sec. \n',t_elapsed_footprints);

clear hazard_intensity % free up memory

% fill some further fields
for event_i=1:n_events
    i1=1+(event_i-1)*(n_prob_events+1);
    i2=  (event_i  )*(n_prob_events+1);
    hazard_info.event_ID(i1:i2)     = event_i*100+(1:n_prob_events+1)-1;
    hazard_info.orig_event_flag(i1) = 1;
    hazard_info.yyyy(i1:i2)         = hazard_yyyy(event_i);
    hazard_info.mm(i1:i2)           = hazard_mm(event_i);
    hazard_info.dd(i1:i2)           = hazard_dd(event_i);
end % event_i

% first complete hazard_info to contain all except intensity
hazard_info.event_count      = n_events_prob;
hazard_info.frequency        = (hazard_info.event_ID*0+1)/(hazard_info.orig_years*(1+n_prob_events));
hazard_info.orig_event_flag  = logical(hazard_info.orig_event_flag); % convert to logical for indexing
hazard_info.t_elapsed_footprints = t_elapsed_footprints; % later to be commented out, for timing
hazard_info.date             = datestr(now);

% Create a yearset to be used in climada_EDS2YDS
if create_yearset
    
    % the beginner does not need to understand whats happening here ;-)
    % see climada_EDS2YDS
    
    year_i=1; % init
    active_year=hazard_info.yyyy(1); % first year
    event_index=[];event_count=0; % init
    
    fprintf('> generating yearset for %i events: ',n_events_prob);
    climada_progress2stdout    % init, see terminate below
    t0       = clock;
    for event_i=1:n_events_prob
        if hazard_info.yyyy(event_i)==active_year
            if hazard_info.orig_event_flag(event_i) % same year, add if original track
                event_count=event_count+1;
                event_index=[event_index event_i];
            end
        else
            % new year, save last year
            hazard_info.orig_yearset(year_i).yyyy=active_year;
            hazard_info.orig_yearset(year_i).event_count=event_count;
            hazard_info.orig_yearset(year_i).event_index=event_index;
            year_i=year_i+1;event_index=[];event_count=0; % re-init
            % reset for next year
            active_year=hazard_info.yyyy(event_i);
            if hazard_info.orig_event_flag(event_i)
                % same year, add if original track
                event_count=1;
                event_index=event_i;
            end
        end
        climada_progress2stdout(event_i,n_events_prob,10,'events'); % update
    end % event_i
    climada_progress2stdout(0) % terminate
    
    % save last year
    hazard_info.orig_yearset(year_i).yyyy=active_year;
    hazard_info.orig_yearset(year_i).event_count=event_count;
    hazard_info.orig_yearset(year_i).event_index=event_index;
    
    t_elapsed_yearset = etime(clock,t0);
    fprintf('done, took %3.2f sec. \n',t_elapsed_yearset);
    
end % create_yearset

% now loop over all countries and write the probabilistic hazard sets
for country_i=n_countries:-1:1 % backard to have last country last in memory
    hazard              = hazard_info; % copy, then restrict:
    hazard              = rmfield(hazard,'shape_i');
    hazard              = rmfield(hazard,'shape');
    hazard.lon          = hazard.lon(country.pos{country_i}.pos);
    hazard.lat          = hazard.lat(country.pos{country_i}.pos);
    hazard.area_km2     = hazard.area_km2(country.pos{country_i}.pos);
    hazard.centroid_ID  = hazard.centroid_ID(country.pos{country_i}.pos);
    hazard.intensity    = hazards{country_i}.intensity;
    if ~strcmpi(climada_global.save_file_version,'-v7') && add_fraction % gets too big in Octave
        hazard.fraction = spones(hazard.intensity); % fraction 100%
    end
    hazard.filename=country.hazard_filename{country_i};
    fprintf('> saving %s\n',hazard.filename);
    save(hazard.filename,'hazard',climada_global.save_file_version);
end % country_i

if check_plot,climada_hazard_plot(hazard,0);end

if nargout>0,hazard_info=hazard;end % to return

end % wisc_hazard_set_prob