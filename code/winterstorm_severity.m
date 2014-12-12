function severity=winterstorm_severity(hazard,check_plot)
% climada
% NAME:
%   winterstorm_severity
% PURPOSE:
%   calculate the winterstorm severity, a measure for it's 'strength'
%
%   in essence, calculate sum of v^3 over all centroids of a each event
%
%   see winterstorm_compare_severity (to show severity exceedence curve etc)
% CALLING SEQUENCE:
%   severity=winterstorm_severity(hazard)
% EXAMPLE:
%   severity=winterstorm_severity(hazard)
% INPUTS:
%   hazard: the filename of a hazard event set
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   check_plot: if =1, show check plot of selected centroids
%       =0: no plot (default)
% OUTPUTS:
%   severity: a structure with
%       index= the storm severity index
%       frequency: the event frequency
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141128, initial
%-

severity=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('hazard','var'),hazard=[];end
if ~exist('check_plot','var'),check_plot=0;end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
country_risk_dir=[fileparts(fileparts(which('country_admin1_risk_calc'))) filesep 'data'];

% PARAMETERS
%
% admin0 shape file (a file in the climada module country_risk (see
% https://github.com/davidnbresch/climada_module_country_risk)
admin0_shape_file=[country_risk_dir filesep 'ne_10m_admin_0_countries'...
    filesep 'ne_10m_admin_0_countries.shp'];
%
% define the countries we would like to calculate the severity for
country_list={'Belgium','Switzerland','Germany','Denmark','France','United Kingdom','Austria','Netherlands','Ireland'};
%
% TEST
%hazard_filename=[module_data_dir filesep 'validation' filesep 'Daria_biasMean.mat'];
%load(hazard_filename)

% prompt for hazard if not given
if isempty(hazard) % local GUI
    hazard_filename=[module_data_dir filesep 'validation' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard_filename, 'Select event set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_filename=fullfile(pathname,filename);
    end
    load(hazard_filename); % load hazard
end

% get the admin0 boundaries
if exist(admin0_shape_file,'file')
    %fprintf('reading admin0 boundaries\n');
    shapes=climada_shaperead(admin0_shape_file); % read the admin0 shapes
    n_shapes=length(shapes);
else
    fprintf('ERROR: admin0 shape file %s not found, aborted\n',admin0_shape_file);
    fprintf('download www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip\n');
    fprintf('and store in data dir of the climada country_risk module (install module first)\n');
    return
end

fprintf('selecting centroids ');
valid_centroids=hazard.lon*0; % init
if climada_global.waitbar,h = waitbar(0,sprintf('selecting from %i countries',n_shapes),'name','Severity: centroid selection');end
for shape_i=1:n_shapes
    if strmatch(shapes(shape_i).NAME,country_list)
        country_centroids=inpolygon(hazard.lon,hazard.lat,shapes(shape_i).X,shapes(shape_i).Y);
        valid_centroids=valid_centroids+country_centroids;
        if climada_global.waitbar,waitbar(shape_i/n_shapes,h);end
    end
end % shape_i
if climada_global.waitbar,close(h);end % dispose waitbar

valid_centroids(find(valid_centroids))=1;
valid_centroids=logical(valid_centroids); % turn to index
n_valid_centroids=sum(valid_centroids);

fprintf('%i (of %i) ',n_valid_centroids,length(valid_centroids));

if check_plot
    climada_plot_world_borders;hold on
    plot(hazard.lon,hazard.lat,'xr');
    plot(hazard.lon(valid_centroids),hazard.lat(valid_centroids),'og');
end

severity.filename=hazard.filename;
severity.frequency=hazard.frequency;
severity.index=severity.frequency*0; % init

msgstr=sprintf('processing %i events at %i centroids',hazard.event_count,n_valid_centroids);
fprintf('%s\n',msgstr);
if climada_global.waitbar,h = waitbar(0,msgstr,'name','Severity: index calculation');end

for event_i=1:hazard.event_count
    event_intensity=full(hazard.intensity(event_i,valid_centroids));
    severity.index(event_i)=sum((event_intensity.^3))/n_valid_centroids;
    severity.index(event_i)=severity.index(event_i)^(1/3);
    
    if climada_global.waitbar && mod(event_i,100)==0
        waitbar(event_i/hazard.event_count,h,msgstr);
    end
    
end % event_i
if climada_global.waitbar,close(h);end % dispose waitbar

return
