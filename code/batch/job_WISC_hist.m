% batch job for cluster: bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC_hist
% MODULE:
%   storm_europe
% NAME:
%   job_WISC_hist
% PURPOSE:
%   run WISC historic WS hazard set generation, in essence a caller
%   for wisc_hazard_set on the cluster, then run some statistics, too.
%
%   See PARAMETERS below before running
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster: scp -r Documents/_GIT/climada_modules/storm_europe/code/job_WISC.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   run on cluster:      bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC_hist
%
%   copy all data to cluster:  
%   scp -r Documents/_GIT/climada_data/WISC/C3S_WISC_FOOTPRINT_NETCDF_0100 dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/WISC/.
%   next line, if it exists, othewrwise consider running WISC_era20c_eur_WS_hist first alone locally before to generate it
%   scp -r Documents/_GIT/climada_data/results/WISC_hazard_plus.mat dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/results/.
%
%   copy results back local:   
%   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/WISC_*hist*.mat Documents/_GIT/climada_data/hazards/.
%   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/results/WISC_*hist*.mat Documents/_GIT/climada_data/results/.
%   copy results back polybox: 
%   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/WISC_*hist*.mat /Users/bresch/polybox/WISC/hazards/.
%   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/results/WISC_*hist*.mat /Users/bresch/polybox/WISC/results/.
%
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC_hist
% EXAMPLE:
%   bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC_hist
%
%   run_on_desktop=1; % to test the job on a desktop
%   job_WISC_hist
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   run_on_desktop: if you set =1 before calling job_WISC_hist (all in
%       MATLAB command window), it sets the number of parallel pool workers to
%       two, does not delete the pool after execution and does not quit
%       MATLAB and hence allows to TEST a job on a local desktop. This
%       allows to TEST a job without editing it.  
%       Default=0 for cluster.
% OUTPUTS:
%   to scratch disk, see PARAMETERS
% MODIFICATION HISTORY:
% David N. Bresch, dbresch@ethz.ch, 20180105, copy from job_WISC
% David N. Bresch, dbresch@ethz.ch, 20180106, run_on_desktop added
%-

% PARAMETERS
%
cluster_climada_root_dir='/cluster/home/dbresch/climada'; % to make sure the cluster finds climada
%
N_pool_workers=24; % for parpool on cluster
%
% define where the WISC windstorm footprints are stored
wisc_dir    = '/cluster/work/climate/dbresch/climada_data/WISC/C3S_WISC_FOOTPRINT_NETCDF_0100/';
%
% pro memoria the key directories on the cluster
%scratch_dir = '/cluster/scratch/dbresch/climada_data/hazards';
%hazards_dir = '/cluster/work/climate/dbresch/climada_data/hazards'; % same as climada_global.hazards_dir if local startupp.m defines it


% some admin to start with
if ~exist('run_on_desktop','var'),run_on_desktop=[];end
if isempty(run_on_desktop),run_on_desktop=0;end % default=0, =1 to run job on mac
if run_on_desktop % for parpool on desktop
    N_pool_workers=2;
    pool_active=gcp('nocreate');
    if isempty(pool_active),pool=parpool(N_pool_workers);end
else
    cd(cluster_climada_root_dir)
    pool=parpool(N_pool_workers);
end
startup % climada_global exists afterwards
fprintf('executing in %s\n',pwd) % just to check where the job is running from
climada_global.parfor=1; % for parpool
t0=clock;


% generate both the era20c and eraint probabilistic hazard sets
% -------------------------------------------------------------
% WARNING: one shall first create the addition of the coastal buffer once,
% then stored in WISC_hazard_plus.mat, see add_on_land in wisc_hazard_set
% code. Otherwise parfor might lead to troubles, as the file
% WISC_hazard_plus.mat might be generated twice (and possibly written to at
% the same time).
%

climada_global_hazards_dir=climada_global.hazards_dir; % parfor prefers this
parfor loop_i=1:2 % a bit a joke, but at least parallel
    
    if loop_i==1
        wisc_file=[wisc_dir filesep 'fp_era20c_*.nc'];
        hazard_set_file=[climada_global_hazards_dir filesep 'WISC_era20c_eur_WS_hist.mat'];
    else
        wisc_file=[wisc_dir filesep 'fp_eraint_*.nc'];
        hazard_set_file=[climada_global_hazards_dir filesep 'WISC_eraint_eur_WS_hist.mat'];
    end
    
    if ~exist(hazard_set_file,'file')
        wisc_hazard_set(wisc_file,0,hazard_set_file);
    else
        fprintf('%s already exists\n',hazard_set_file);
    end
    
end % loop_i


% merge the two hazard sets
% -------------------------

hazard_set_file=[climada_global.hazards_dir filesep 'WISC_eur_WS_hist.mat'];
if ~exist(hazard_set_file,'file')
    load([climada_global.hazards_dir filesep 'WISC_eraint_eur_WS_hist.mat']); % load second hazard first, since larger
    hazard_eraint=hazard;clear hazard
    load([climada_global.hazards_dir filesep 'WISC_era20c_eur_WS_hist.mat']); % load first hazard
    climada_hazard_merge(hazard,hazard_eraint,'events',hazard_set_file);
else
    fprintf('%s already exists\n',hazard_set_file)
end


% calculate statistics
% --------------------

clear hazard % free up memory
global hazard % define as global to allow for functions to work on without passing on the whole data
load([climada_global.hazards_dir filesep 'WISC_eur_WS_hist.mat']); % load hazard

% calculate storm severity index (SSI)
ssi_filename=[climada_global.results_dir filesep 'WISC_eur_WS_hist_ssi.mat'];
[ssi.ssi,ssi.ssi_sorted,ssi.xs_freq,ssi.ssi_orig,ssi.ssi_sorted_orig,ssi.xs_freq_orig]=climada_hazard_ssi('global',0);
fprintf('> storing ssi to %s\n',ssi_filename);
save(ssi_filename,'ssi',climada_global.save_file_version)

% calculate hazard return period maps
return_periods=[5 10 25 50 100 250 500 1000]; % to be consistent

stats_filename=[climada_global.results_dir filesep 'WISC_eur_WS_hist_stats.mat'];
climada_hazard_stats('global',return_periods,-10); % operate on global variable hazard
stats_hist.stats       = hazard.stats;
stats_hist.lon         = hazard.lon;
stats_hist.lat         = hazard.lat;
stats_hist.peril_ID    = hazard.peril_ID;
stats_hist.units       = hazard.units;
stats_hist.event_ID    = hazard.event_ID;
stats_hist.centroid_ID = hazard.centroid_ID;
		
fprintf('> storing historic stats to %s\n',stats_filename);
save(stats_filename,'stats_hist',climada_global.save_file_version)

climada_hazard_stats('global',return_periods,  0); % operate on global variable hazard
stats.stats       = hazard.stats;
stats.lon         = hazard.lon;
stats.lat         = hazard.lat;
stats.peril_ID    = hazard.peril_ID;
stats.units       = hazard.units;
stats.event_ID    = hazard.event_ID;
stats.centroid_ID = hazard.centroid_ID;
fprintf('> storing historic and probabilistic stats to %s\n',stats_filename);
save(stats_filename,'stats','stats_hist',climada_global.save_file_version)


% visualize the statistics
% ------------------------
clear hazard % free up memory

% plot storm severity index (SSI)
load([climada_global.results_dir filesep 'WISC_eur_WS_hist_ssi.mat']);
fig_ssi=figure('Name','WISC hazard ssi','Visible','off');
plot(ssi.xs_freq,ssi.ssi_sorted);
xlabel('xs frequency (years)');ylabel('SSI (arbitrary units)');set(gcf,'Color',[1 1 1]);hold on
plot(ssi.xs_freq_orig,ssi.ssi_sorted_orig,'xr');
legend({'all events','historic only'});
saveas(fig_ssi,[climada_global.results_dir filesep 'WISC_eur_WS_hist_ssi'],'png');delete(fig_ssi); % save and delete figure

% plot hazard return period maps
load([climada_global.results_dir filesep 'WISC_eur_WS_hist_stats.mat']);
hazard.stats       = stats_hist.stats;
hazard.lon         = stats_hist.lon;
hazard.lat         = stats_hist.lat;
hazard.peril_ID    = stats_hist.peril_ID;
hazard.units       = stats_hist.units;
hazard.event_ID    = stats_hist.event_ID;
hazard.centroid_ID = stats_hist.centroid_ID;
hazard.intensity   = []; % dummy
hazard.frequency   = []; % dummy
fig_stats=figure('Name','WISC hazard stats','Visible','off');
climada_hazard_stats(hazard,hazard.stats.return_period,[-1,4,2]); % historic, 4x2 (hor x vert)
saveas(fig_stats,[climada_global.results_dir filesep 'WISC_eur_WS_hist_stats_hist'],'png');delete(fig_stats); % save and delete figure

load([climada_global.results_dir filesep 'WISC_eur_WS_hist_stats.mat']);
hazard.stats       = stats.stats;
hazard.lon         = stats.lon;
hazard.lat         = stats.lat;
hazard.peril_ID    = stats.peril_ID;
hazard.units       = stats.units;
hazard.event_ID    = stats.event_ID;
hazard.centroid_ID = stats.centroid_ID;
hazard.intensity   = []; % dummy
hazard.frequency   = []; % dummy
fig_stats=figure('Name','WISC hazard stats','Visible','off');
climada_hazard_stats(hazard,hazard.stats.return_period,[ 1,4,2]); % probailistic, 4x2 (hor x vert)
saveas(fig_stats,[climada_global.results_dir filesep 'WISC_eur_WS_hist_stats'],'png');delete(fig_stats); % save and delete figure


% some admin to close
fprintf('total job time %f seconds\n',etime(clock, t0));
if ~run_on_desktop,delete(pool);exit;end % no need to delete the pool on mac, the cluster appreciates exit