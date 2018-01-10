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
%   previous call: <note the most usual previous call here>
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
%       for details), just for the combined WISC hazards. 
%   combined_ssi: the combined Storm Severity Index (SSI)
%   country: the country information, to still know which part belongs to
%       which country, a struct with fields:
%   shape_i: the shape number (of ID) for each centroid, i.e. to which
%       country shape each centroid belongs to. Assigned in
%       wisc_hazard_set_prob, see hazard_info.filename in PARAMETERS
%   shape: a struct with fields:
%       ISO3{i}: the ISO3 code for country i
%       shape_i{i}: the ID in shape_i for country i       
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

load(hazard_info.filename)
n_centroids=length(hazard_info.lon);

if ~isempty(return_periods) % init
    n_return_periods=length(return_periods);
    combined_stats.intensity=spalloc(n_return_periods,n_centroids,ceil(n_return_periods*n_centroids*matrix_density));
    combined_stats.shape_i       = hazard_info.shape_i;
    combined_stats.shape         = hazard_info.shape;
end

country_i=1;
for file_i=1:length(wisc_files)
    clear hazard;hazard=[];
    hazard_filename=[wisc_dir filesep wisc_files(file_i).name];
    load(hazard_filename); % loads hazard
    
    if ~isempty(hazard)
        
        [~,fN]=fileparts(hazard_filename);
        fprintf('> adding %s\n',fN);
        
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
                combined_stats.lon=hazard_info.lon;
                combined_stats.lat=hazard_info.lat;
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
            if isempty(combined_ssi)
                combined_ssi      = ssi;
            else
                combined_ssi      = combined_ssi + ssi;
            end
        end % calc_ssi>0
        
    end % ~isempty(hazard)
    
end % file_i

% plotting
% --------
if check_plot
    if ~isempty(combined_stats)
        climada_hazard_stats(combined_stats,return_periods,1);
    end
    
    if ~isempty(combined_ssi)
        
        [ssi_sorted,xs_freq]=climada_damage_exceedence(combined_ssi,hazard.frequency,hazard.event_ID,1);
        
        if isfield(hazard,'orig_event_flag')
            hazard_orig_event_flag=logical(hazard.orig_event_flag); % to be sure
            ssi_orig=combined_ssi(hazard_orig_event_flag);
            frequency_orig=hazard.frequency(hazard_orig_event_flag)*hazard.event_count/hazard.orig_event_count;
            [ssi_sorted_orig,xs_freq_orig]=climada_damage_exceedence(ssi_orig,frequency_orig,[],1);
        end
        
        %plot(return_period,ssi_sorted);xlabel('return period (years)');ylabel('SSI');set(gcf,'Color',[1 1 1])
        plot(xs_freq,ssi_sorted);
        xlabel('xs frequency (years)');ylabel('SSI (arbitrary units)');set(gcf,'Color',[1 1 1]);hold on
        
        if isfield(hazard,'orig_event_flag')
            plot(xs_freq_orig,ssi_sorted_orig,'xr');
            legend({'all events','historic only'});
        end
        
    end
end % check_plot

end % wisc_hazard_stats