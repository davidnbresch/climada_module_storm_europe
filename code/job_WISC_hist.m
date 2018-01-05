% batch job for cluster: bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC_hist
% MODULE:
%   storm_europe
% NAME:
%   job_WISC_hist
% PURPOSE:
%   run WISC probabilistic WS hazard set generation, in essence a caller
%   for wisc_hazard_set on the cluster, then run some statistics, too.
%
%   See PARAMETERS below before running
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r Documents/_GIT/climada_modules/storm_europe/code/job_WISC_hist.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%
%   copy all data to cluster:  scp -r Documents/_GIT/climada_data/WISC/C3S_WISC_FOOTPRINT_NETCDF_0100 dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/WISC/.
%   copy all data to cluster:  scp -r Documents/_GIT/climada_data/hazards/WISC_eur_WS_hist.mat dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/.
%   copy all data to cluster:  scp -r Documents/_GIT/climada_data/results/WISC_hazard_plus.mat dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/results/.
%   run on cluster:            bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC_hist
%
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/WISC_eur_WS*.mat Documents/_GIT/climada_data/hazards/.
%   copy results back polybox: scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/WISC_eur_WS*.mat /Users/bresch/polybox/WISC/hazards/.
%
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC_hist
% EXAMPLE:
%   bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC_hist
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to scratch disk, see PARAMETERS
% MODIFICATION HISTORY:
% David N. Bresch, dbresch@ethz.ch, 20171230, copy from job_isimip04
% David N. Bresch, dbresch@ethz.ch, 20180101, run first, output to /cluster/work/...
% David N. Bresch, dbresch@ethz.ch, 20180102, run again for WISC_eraint_eur_WS
% David N. Bresch, dbresch@ethz.ch, 20180103, run again for merging
% David N. Bresch, dbresch@ethz.ch, 20180105, run again for stats
%-

% PARAMETERS
%
% one should only have to edit this section
cd /cluster/home/dbresch/climada % to make sure the cluster finds climada
%scratch_dir = '/cluster/scratch/dbresch/climada_data/hazards';
%wisc_dir   ='/cluster/work/climate/dbresch/climada_data/WISC/C3S_WISC_FOOTPRINT_NETCDF_0100/';
%hazards_dir='/cluster/work/climate/dbresch/climada_data/hazards'; %same as climada_global.hazards_dir if local startupp.m defines it


startup % climada_global exists afterwards
pwd % just to check where the job is running from
N_pool_workers=24; % for parpool
climada_global.parfor=1; % for parpool


pool=parpool(N_pool_workers);

% generate both the era20c and eraint probabilistic hazard sets
% -------------------------------------------------------------
% WARNING: one shall first create the addition of the coastal buffer once,
% then stored in WISC_hazard_plus.mat, see add_on_land in wisc_hazard_set
% code. Otherwise parfor might lead to troubles, as the file
% WISC_hazard_plus.mat might be generated twice (and possibly written to at
% the same time).
%
% parfor loop_i=1:2 % a bit a joke, but at least parallel
%
%     if loop_i==1
%         hazard_era20c=wisc_hazard_set([wisc_dir filesep 'fp_era20c_*.nc'],0,'WISC_era20c_eur_WS_hist');
%     else
%         hazard_eraint=wisc_hazard_set([wisc_dir filesep 'fp_eraint_*.nc'],0,'WISC_eraint_eur_WS_hist');
%     end
%
% end % loop_i

% merge the two hazard sets
% -------------------------

% load([hazards_dir filesep 'WISC_eraint_eur_WS_hist.mat']); % load first hazard
% hazard_eraint=hazard;clear hazard
% load([hazards_dir filesep 'WISC_era20c_eur_WS_hist.mat']); % load second hazard
%
% climada_hazard_merge(hazard,hazard_eraint,'events',[hazards_dir filesep 'WISC_eur_WS_hist' '.mat']);

% calculate statistics
% --------------------

clear hazard % free up memory
global hazard % define as global to allow for functions to work on without passing on the whole data
load([climada_global.hazards_dir filesep 'WISC_eur_WS_hist.mat']); % load hazard

% calculate storm severity index (SSI)
ssi_filename=[climada_global.hazards_dir filesep 'WISC_eur_WS_hist_ssi.mat'];
tic
[ssi.ssi,ssi.ssi_sorted,ssi.xs_freq,ssi.ssi_orig,ssi.ssi_sorted_orig,ssi.xs_freq_orig]=climada_hazard_ssi('global',0);
toc
fprintf('> storing ssi to %s\n',ssi_filename);
save(ssi_filename,'ssi',climada_global.save_file_version)

% calculate hazard return period maps
stats_filename=[climada_global.hazards_dir filesep 'WISC_eur_WS_hist_stats.mat'];
tic
climada_hazard_stats('global',[10 25 50 100 250 500 1000],-10); % operate on global variable hazard
toc
stats_hist = hazard.stats;
fprintf('> storing historic stats to %s\n',stats_filename);
save(stats_filename,'stats_hist',climada_global.save_file_version)

tic
climada_hazard_stats('global',[10 25 50 100 250 500 1000],-10); % operate on global variable hazard
toc
stats      = hazard.stats;
fprintf('> storing historic and probabilistic stats to %s\n',stats_filename);
save(stats_filename,'stats','stats_hist',climada_global.save_file_version)

delete(pool)

exit % the cluster appreciates this, gives back memory etc.