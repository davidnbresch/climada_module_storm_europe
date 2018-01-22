function [combined_stats,combined_ssi,country]=wisc_hazard_stats(return_periods,calc_ssi,check_plot)
% climada WISC wisc hazard stats
% MODULE:
%   storm_europe
% NAME:
%   wisc_hazard_stats
% PURPOSE:
%   compile statistics of wisc hazard sets by compositing country results
%   to whole Europe
%
%   NOTE: this code listens to climada_global.parfor=1 for substantial speedup
%
%   previous call: wisc_hazard_set_prob
%   next call: <note the most usual next function call here>
% CALLING SEQUENCE:
%   [combined_stats,combined_ssi,country]=wisc_hazard_stats(return_periods,calc_ssi,check_plot)
% EXAMPLE:
%   wisc_hazard_set_prob % be careful, this takes time
%   [combined_stats,combined_ssi]=wisc_hazard_stats([10 25 50 100 250 1000],1,1)
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   return_periods: the return periods for return period maps. No maps if
%       empty
%   calc_ssi: whether we calc Storm Severity Index (SSI), =1 (default), or
%       not (=0)
%   check_plot: if =1, plot the return period maps and SSI distribution, default =0 (no plot)
% OUTPUTS:
%   combined_stats: the same output as from climada_hazard_stats (see there
%       for details), just for the combined WISC hazards. Key fields are:
%       historic: =1 if for historic events only, see check_plot=-1 or -10
%       return_period(i): return period i
%       intensity(i,j): intensity for return period i at centroid j
%       One can easily plot single return period maps using
%       climada_hazard_plot(combined_stats,1) where 1 for the first return
%       period etc. and then label e.g. as
%       title(sprintf('%i yr hazard map',combined_stats.return_period(1)))
%   combined_ssi: the Storm Severity Index (SSI), a struct with fields:
%       ssi: the SSI values (for each event)
%       ssi_sorted: sorted ssi
%       xs_freq: the exceedance frequency, to plot(xs_freq,ssi_sorted)
%       *_hist: the same, but for historic events only
%   country: the country information, to still know which part belongs to
%       which country, a struct with fields:
%       shape_i: the shape number (of ID) for each centroid, i.e. to which
%           country shape each centroid belongs to. Assigned in
%           wisc_hazard_set_prob, see hazard_info.filename in PARAMETERS
%       shape: a struct with fields:
%           ISO3{i}: the ISO3 code for country i
%           shape_i{i}: the ID in shape_i for country i
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20180110, initial
%-

combined_stats=[];combined_ssi=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('return_periods','var'),return_periods=[];end
if ~exist('calc_ssi','var'),      calc_ssi=1;end
if ~exist('check_plot','var'),    check_plot=0;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
hazard_info.filename=[climada_global.data_dir filesep 'WISC' filesep 'WISC_hazard_info.mat'];
%
matrix_density=.4; % the average density of the hazard (=#nozero/#all)
%
WISC_hazard_files=[climada_global.hazards_dir filesep 'WISC_*_eur_WS.mat'];


wisc_dir=fileparts(WISC_hazard_files);
wisc_files=dir(WISC_hazard_files);

fprintf('< loading %s\n',hazard_info.filename);
load(hazard_info.filename)
n_centroids=length(hazard_info.lon);

if ~isempty(return_periods) % init
    n_return_periods=length(return_periods);
    combined_stats.intensity=spalloc(n_return_periods,n_centroids,ceil(n_return_periods*n_centroids*matrix_density));
    combined_stats.shape_i       = hazard_info.shape_i;
    combined_stats.shape         = hazard_info.shape;
end

country_i=1;ssi_combined=[];
for file_i=1:length(wisc_files)
%for file_i=1:3
    clear hazard;hazard=[];
    hazard_filename=[wisc_dir filesep wisc_files(file_i).name];
    load(hazard_filename); % loads hazard
    
    if ~isempty(hazard)
        
        [~,fN]=fileparts(hazard_filename);
        fprintf('< adding %s\n',fN);
        
        if ~isempty(return_periods)
            stats=climada_hazard_stats(hazard,return_periods,0);
            
            country.ISO3{country_i}=strrep(strrep(fN,'WISC_',''),'_eur_WS','');
            iii=[]; % figure idex of centroids for each country (for speedup)
            for si=1:length(hazard_info.shape.shape_i)
                % find country.ISO3 in hazard_info.shape.ISO3
                if strcmpi(hazard_info.shape.ISO3{si},country.ISO3{country_i}),iii=si;end
            end %si
            country.pos{country_i}.pos=find(hazard_info.shape_i==hazard_info.shape.shape_i(iii)); % returns index values, saves space
            country_i=country_i+1;
            % e.g.: hazard.shape_i(country.pos{country_i}.pos)) % all centroids within country i
            
            if file_i==1 % init
                combined_stats.lon           = hazard_info.lon;
                combined_stats.lat           = hazard_info.lat;
                combined_stats.centroid_ID   = 1:n_centroids;
                combined_stats.historic      = stats.historic;
                combined_stats.return_period = stats.return_period;
                combined_stats.peril_ID      = stats.peril_ID;
                combined_stats.units         = stats.units;
                combined_stats.frequency     = stats.frequency;
                combined_stats.event_ID      = stats.event_ID;
                combined_stats.shape_i       = hazard_info.shape_i;
                combined_stats.shape         = hazard_info.shape;
            end
            combined_stats.intensity(:,country.pos{country_i-1}.pos)=stats.intensity;
            
        end % ~isempty(return_periods)
        
        if calc_ssi>0
            ssi=climada_hazard_ssi(hazard,0);
            if isempty(ssi_combined)
                ssi_combined      = ssi;
            else
                ssi_combined      = ssi_combined + ssi;
            end
        end % calc_ssi>0
        
    end % ~isempty(hazard)
    
end % file_i

% process ssi
% -----------

[ssi_sorted,xs_freq]=climada_damage_exceedence(ssi_combined,hazard.frequency,hazard.event_ID,1);
clear ssi
combined_ssi.ssi=ssi_combined;
combined_ssi.ssi_sorted=ssi_sorted;
combined_ssi.xs_freq=xs_freq;

if isfield(hazard,'orig_event_flag') % orig_event_flag just on last hazard
    hazard_orig_event_flag=logical(hazard.orig_event_flag); % to be sure
    ssi_orig=ssi_combined(hazard_orig_event_flag);
    frequency_orig=hazard.frequency(hazard_orig_event_flag)*hazard.event_count/hazard.orig_event_count;
    [ssi_sorted_orig,xs_freq_orig]=climada_damage_exceedence(ssi_orig,frequency_orig,[],1);
    combined_ssi.ssi_orig=ssi_orig;
    combined_ssi.ssi_sorted_orig=ssi_sorted_orig;
    combined_ssi.xs_freq_orig=xs_freq_orig;
end

% plotting
% --------
if check_plot
    if ~isempty(combined_stats)
        figure('Name','return period maps');
        climada_hazard_stats(combined_stats,return_periods,1);
    end
    
    if ~isempty(ssi_combined)
        
        figure('Name','SSI');
        plot(xs_freq,ssi_sorted);
        xlabel('xs frequency (years)');ylabel('SSI (arbitrary units)');set(gcf,'Color',[1 1 1]);hold on
        
        if isfield(hazard,'orig_event_flag')
            plot(xs_freq_orig,ssi_sorted_orig,'xr');legend({'all events','historic only'});
        end
        
    end
end % check_plot

end % wisc_hazard_stats