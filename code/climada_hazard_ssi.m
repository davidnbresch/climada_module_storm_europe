function [ssi,sorted_ssi,exceedence_freq]=climada_hazard_ssi(hazard,check_plot,windspeed_threshold_ms)
% climada hazard winter storm ws ssi
% MODULE:
%   storm_europe
% NAME:
%   climada_hazard_ssi
% PURPOSE:
%   Calculate Storm Severity Index (SSI) according to Lamb (1990)
%
%   SSI = v [m/s] ^ 3 * duration [h] * area [km^2 or m^2], here with duration =1
%   summed up over one event for all centroids over land where v >= windspeed_threshold_ms
%
%   previous call: generate a winter storm hazard event set, e.g. see wisc_hazard_set
%   next call: muse about the output
% CALLING SEQUENCE:
%   ssi=climada_hazard_ssi(hazard,check_plot,windspeed_threshold_ms)
% EXAMPLE:
%   hazard=wisc_hazard_set('test'); % create a windstorm hazard
%   ssi=climada_hazard_ssi(hazard,1);
% INPUTS:
%   hazard: a climada hazard set, especially useful for winter storm. Should
%       contain the field hazard.area_km2, the area per centroid (otherwise, a
%       uniform area of 1 is presumed). Should also contab the field
%       hazard.on_land, as otherwise all centroids are used.
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   check_plot: if =1, plot the SSI distribution, default =0 (no plot)
%   windspeed_threshold_ms: the windspeed threshold in m/s, default=22 (Lamb)
% OUTPUTS:
%   ssi(i): the ssi index for each event i
%   sorted_ssi(j): same as ssi, sorted descendingly
%   exceedence_freq(j): xs frequency for sorted_ssi(j). 
%       Useful e.g. to plot(exceedence_freq,sorted_ssi);
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20171229, initial
%-

ssi=[]; % init output

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('hazard','var'),hazard=[];end
if ~exist('check_plot','var'),check_plot=0;end
if ~exist('windspeed_threshold_ms','var'),windspeed_threshold_ms=22;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%

hazard = climada_hazard_load(hazard); % prompt for hazard_set if not given
if isempty(hazard),return;end
hazard = climada_hazard2octave(hazard); % Octave compatibility for -v7.3 mat-files

if ~isfield(hazard,'area_km2')
    fprintf('WARNING: no field hazard.area_km2, unit area presumed\n');
    hazard.area_km2=hazard.lon*0+1; % unit area presumed
end

if isfield(hazard,'on_land')
    hazard.area_km2=hazard.area_km2.*hazard.on_land; % deal with land points only
else
    fprintf('WARNING: no field hazard.on_land, all centroids used\n');
end
 
n_events=size(hazard.intensity,1);
n_centroids=size(hazard.intensity,2);
fprintf('> calculating storm severity index for %i events (at %i centroids): ',n_events,n_centroids);

climada_progress2stdout    % init, see terminate below
t0 = clock;

ssi=zeros(1,n_events); % init
for event_i=1:n_events
    [~,cols,intensity] = find(hazard.intensity(event_i,:));
    intensity(intensity<windspeed_threshold_ms)=0;
    ssi(event_i) = intensity.^3 * hazard.area_km2(cols)'*1e6; % v^3*m^2=(m/s)^3*m^2=m^5/s^5
    climada_progress2stdout(event_i,n_events,10,'events'); % update
end % event_i

% % slower
% for event_i=1:n_events
%     [~,cols] = find(hazard.intensity(event_i,:)>=windspeed_threshold_ms);
%     ssi(event_i) = hazard.intensity(event_i,cols).^3 * hazard.area_km2(cols)'*1e6; % v^3*m^2=(m/s)^3*m^2=m^5/s^5
%     climada_progress2stdout(event_i,n_events,10,'events'); % update
% end % event_i

% % slower
% hazard.intensity(hazard.intensity<windspeed_threshold_ms)=0;
% hazard.intensity=hazard.intensity.^3; % cube
% ssi=hazard.area_km2*hazard.intensity'*1e6; % v^3*m^2=(m/s)^3*m^2=m^5/s^5

climada_progress2stdout(0) % terminate

t_elapsed = etime(clock,t0);
fprintf('done, took %3.2f sec. \n',t_elapsed);

% scale ssi
ssi=ssi*1e-12; % arbitrary scaling
% 8e-9 from old approach, close to Lamb

[sorted_ssi,exceedence_freq]=climada_damage_exceedence(ssi,hazard.frequency,hazard.event_ID,1);

if check_plot
    %return_period   = 1./exceedence_freq;
    %plot(return_period,sorted_ssi);xlabel('return period (years)');ylabel('SSI');set(gcf,'Color',[1 1 1])
    plot(exceedence_freq,sorted_ssi);
    xlabel('xs frequency (years)');ylabel('SSI (arbitrary units)');set(gcf,'Color',[1 1 1]);hold on
    
    if isfield(hazard,'orig_event_flag')
        ssi_orig=ssi(logical(hazard.orig_event_flag));
        %ssi_orig=ssi(hazard.orig_event_flag);
        frequency_orig=hazard.frequency(hazard.orig_event_flag)*hazard.event_count/hazard.orig_event_count;
        [sorted_ssi_orig,exceedence_freq_orig]=climada_damage_exceedence(ssi_orig,frequency_orig,[],1);
        plot(exceedence_freq_orig,sorted_ssi_orig,'xr');
        legend({'all events','historic only'});
    end

%     % produce gamma fit
%     % -----------------
% 
%     p_alpha=0.0014; % alpha-parameter of Gamma-distribution (effect: mainly in SSI-direction)
%     all_freq_sum=sum(hazard.frequency);
%     p_gamma=hazard.event_count/(max(hazard.yyyy)-min(hazard.yyyy)+1);
%     
%     ssi_truncation=min(ssi)/2;
%     max_axis_ssi=max(ssi)*1.5;
%     fprintf('Gamma-fit parameters alpha=%f gamma=%f SSI truncation=%f\n',...
%         p_alpha,p_gamma,ssi_truncation);
%     ssi_axis_step=max_axis_ssi/100; % resolution of SSI axis
%     ssi_axis=(ssi_truncation:ssi_axis_step:max_axis_ssi);
%     
%     gamma_fit=p_alpha^p_gamma/gamma(p_gamma)*(ssi_axis-ssi_truncation).^(p_gamma-1).*...
%         exp(-p_alpha*(ssi_axis-ssi_truncation));
%     gamma_fit=gamma_fit*ssi_axis_step; % not done all operations in one to stay close to EXCEL sheet 
%     gamma_fit=cumsum(gamma_fit);
%     gamma_fit=1-gamma_fit;
%     gamma_fit=all_freq_sum*gamma_fit;
%     
%     % make the plot
%     % -------------
%     plot(gamma_fit,ssi_axis,'.k','MarkerSize',1); % overlay gamma fit
%     legend({'SSI','Gamma fit'});
    
end % check_plot

end % climada_hazard_ssi