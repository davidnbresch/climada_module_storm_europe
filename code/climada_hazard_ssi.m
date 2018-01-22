function [ssi,ssi_sorted,xs_freq,ssi_orig,ssi_sorted_orig,xs_freq_orig]=climada_hazard_ssi(hazard,check_plot,windspeed_threshold_ms)
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
%   SSI will be calculated both for all and historic events
%
%   NOTE: this code listens to climada_global.parfor for substantial speedup
%
%   previous call: generate a winter storm hazard event set, e.g. see wisc_hazard_set
%   next call: muse about the output
% CALLING SEQUENCE:
%   ssi=climada_hazard_ssi(hazard,check_plot,windspeed_threshold_ms)
% EXAMPLE:
%   hazard=wisc_hazard_set('test'); % create a windstorm hazard
%   ssi=climada_hazard_ssi(hazard,1);
%
%   global hazard % define hazard as global variable
%   load('WISC_eur_WS'); % better use load than climada_hazard_load
%   climada_hazard_ssi('global',1); % operate on global variable hazard
% INPUTS:
%   hazard: a climada hazard set, especially useful for winter storm. Should
%       contain the field hazard.area_km2, the area per centroid (otherwise, a
%       uniform area of 1 is presumed). Should also contab the field
%       hazard.on_land, as otherwise all centroids are used.
%       ='global': (DISABLED) define hazard as a global variable before calling, in
%       which case its workspace is shared and prevents from passing on a
%       huge variable, in which case MATLAB makes a copy (and hence starts
%       swapping for huge hazard sets...). In this case, do NOT specify any
%       ouput (as otherwise, there might be troubles in return).
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   check_plot: if =1, plot the SSI distribution, default =0 (no plot)
%       If climada_hazard_ssi is called with no output arguments and
%       check:plot is not defined, it set set =1 by default.
%   windspeed_threshold_ms: the windspeed threshold in m/s, default=22 (Lamb)
% OUTPUTS:
%   ssi(i): the ssi index for each event i
%   ssi_sorted(j): same as ssi, sorted descendingly
%   xs_freq(j): xs frequency for ssi_sorted(j).
%       Useful e.g. to plot(xs_freq,ssi_sorted);
%   *_orig: same as above, but for original events only (historic ones)
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20171229, initial
% David N. Bresch, david.bresch@gmail.com, 20171230, ssi_orig added
% David N. Bresch, david.bresch@gmail.com, 20180105, hazard='global' added, parfor
% David N. Bresch, david.bresch@gmail.com, 20180109, no change to hazard in order to allow passing as reference
%-

ssi=[];ssi_sorted=[];xs_freq=[];ssi_orig=[];ssi_sorted_orig=[];xs_freq_orig=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('hazard','var'),hazard=[];end
if ~exist('check_plot','var'),check_plot=[];end
if ~exist('windspeed_threshold_ms','var'),windspeed_threshold_ms=22;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%

if isempty(check_plot)
    if nargout==0 
        check_plot=1;
    else
        check_plot=0;
    end
end
% if strcmpi(hazard,'global')
%     % see PROGRAMMERS APOLOGY in header
%     clear hazard
%     fprintf('working on global hazard\n')
%     global hazard % the parser does not like this, we know ;-)
% else
%     hazard=climada_hazard_load(hazard);
% end

if isempty(hazard),return;end
hazard = climada_hazard2octave(hazard); % Octave compatibility for -v7.3 mat-files

if ~isfield(hazard,'area_km2')
    fprintf('WARNING: no field hazard.area_km2, unit area presumed\n');
    hazard_area_km2=hazard.lon*0+1; % unit area presumed
else
    hazard_area_km2=hazard.area_km2; % make direct (speedup)
end

if isfield(hazard,'on_land')
    hazard_area_km2=hazard_area_km2.*hazard.on_land; % deal with land points only
else
    fprintf('WARNING: no field hazard.on_land, all centroids used\n');
end

n_events=size(hazard.intensity,1);
n_centroids=size(hazard.intensity,2);

ssi=zeros(1,n_events); % init

t0 = clock;

if climada_global.parfor
    
    fprintf('> calculating storm severity index for %i events (at %i centroids, parfor) ',n_events,n_centroids);

    intensity=hazard.intensity; % as parfor does not like structs
    %frequency=hazard.frequency;

    parfor event_i=1:n_events
        [~,cols,intensity_tmp] = find(intensity(event_i,:));
        % next line does not work, as it returns a logical array, not the values in intensity_tmp
        %[~,cols,intensity_tmp]=find(intensity(event_i,:)>=windspeed_threshold_ms); % is wrong
        intensity_tmp(intensity_tmp<windspeed_threshold_ms)=0;
        ssi(event_i) = intensity_tmp.^3 * hazard_area_km2(cols)'*1e6; % v^3*m^2=(m/s)^3*m^2=m^5/s^5
    end % event_i
else
    
    fprintf('> calculating storm severity index for %i events (at %i centroids): ',n_events,n_centroids);

    climada_progress2stdout    % init, see terminate below
    n_event_update=max(10^floor(log10(4123))/10,10);
    
    for event_i=1:n_events
        [~,cols,intensity] = find(hazard.intensity(event_i,:));
        intensity(intensity<windspeed_threshold_ms)=0;
        ssi(event_i) = intensity.^3 * hazard_area_km2(cols)'*1e6; % v^3*m^2=(m/s)^3*m^2=m^5/s^5
        climada_progress2stdout(event_i,n_events,n_event_update,'events'); % update
    end % event_i
    
    % % slower
    % for event_i=1:n_events
    %     [~,cols] = find(hazard.intensity(event_i,:)>=windspeed_threshold_ms);
    %     ssi(event_i) = hazard.intensity(event_i,cols).^3 * hazard_area_km2(cols)'*1e6; % v^3*m^2=(m/s)^3*m^2=m^5/s^5
    %     climada_progress2stdout(event_i,n_events,10,'events'); % update
    % end % event_i
    
    % % slower
    % hazard.intensity(hazard.intensity<windspeed_threshold_ms)=0;
    % hazard.intensity=hazard.intensity.^3; % cube
    % ssi=hazard_area_km2*hazard.intensity'*1e6; % v^3*m^2=(m/s)^3*m^2=m^5/s^5
    
    climada_progress2stdout(0) % terminate
end

t_elapsed = etime(clock,t0);
fprintf('done, took %3.2f sec. \n',t_elapsed);

% scale ssi
ssi=ssi*1e-12; % arbitrary scaling
% 8e-9 from old approach, close to Lamb

if nargout==1 && check_plot==0,return;end % only ssi requested

[ssi_sorted,xs_freq]=climada_damage_exceedence(ssi,hazard.frequency,hazard.event_ID,1);

if isfield(hazard,'orig_event_flag')
    hazard_orig_event_flag=logical(hazard.orig_event_flag); % to be sure
    ssi_orig=ssi(hazard_orig_event_flag);
    frequency_orig=hazard.frequency(hazard_orig_event_flag)*hazard.event_count/hazard.orig_event_count;
    [ssi_sorted_orig,xs_freq_orig]=climada_damage_exceedence(ssi_orig,frequency_orig,[],1);
end

if check_plot
    %return_period   = 1./xs_freq;
    %plot(return_period,ssi_sorted);xlabel('return period (years)');ylabel('SSI');set(gcf,'Color',[1 1 1])
    plot(xs_freq,ssi_sorted);
    xlabel('xs frequency (years)');ylabel('SSI (arbitrary units)');set(gcf,'Color',[1 1 1]);hold on
    
    if isfield(hazard,'orig_event_flag')
        plot(xs_freq_orig,ssi_sorted_orig,'xr');
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