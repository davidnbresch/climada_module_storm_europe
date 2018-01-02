% batch job for cluster: bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC
% MODULE:
%   storm_europe
% NAME:
%   job_WISC
% PURPOSE:
%   run WISC probabilistic WS hazard set generation, in essence a caller
%   for wisc_hazard_set on the cluster
%
%   See PARAMETERS below before running
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r Documents/_GIT/climada_modules/storm_europe/code/job_WISC.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   copy all data to cluster:  scp -r Documents/_GIT/climada_data/WISC/C3S_WISC_FOOTPRINT_NETCDF_0100 dbresch@euler.ethz.ch:/cluster/home/dbresch/climada_data/WISC/.
%   copy all data to cluster:  scp -r Documents/_GIT/climada_data/results/WISC_hazard_plus.mat dbresch@euler.ethz.ch:/cluster/home/dbresch/climada_data/results/.
%   run on cluster:            bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC
%
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/scratch/dbresch/climada_data/hazards/WISC_eur_WS.mat Documents/_GIT/climada_data/hazards/.
%   copy results back polybox: scp -r dbresch@euler.ethz.ch:/cluster/scratch/dbresch/climada_data/hazards/WISC_eur_WS.mat /Users/bresch/polybox/WISC/hazards/.
%
%   in the (very unlikely) case you need to copy the climada code, too:
%       scp -r Documents/_GIT/climada_modules/storm_europe/code dbresch@euler.ethz.ch:/cluster/home/dbresch/climada_modules/storm_europe/.
%       scp -r Documents/_GIT/climada/code dbresch@euler.ethz.ch:/cluster/home/dbresch/climada/.
%
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC
% EXAMPLE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_WISC
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to scratch disk, see PARAMETERS
% MODIFICATION HISTORY:
% David N. Bresch, dbresch@ethz.ch, 20171230, copy from job_isimip04
%-

% PARAMETERS
%
% one should only have to edit this section
cd /cluster/home/dbresch/climada % to make sure the cluster finds climada
%scratch_dir = '/cluster/scratch/dbresch/climada_data/hazards';
wisc_dir   ='/cluster/work/climate/dbresch/climada_data/WISC/C3S_WISC_FOOTPRINT_NETCDF_0100/';
hazards_dir='/cluster/work/climate/dbresch/climada_data/hazards';


startup % climada_global exists afterwards
pwd % just to check where the job is running from
N_pool_workers=24; % for parpool
climada_global.parfor=1; % for parpool


pool=parpool(N_pool_workers);

for loop_i=1:2
    
    if loop_i==1
        hazard_era20c=wisc_hazard_set([wisc_dir filesep 'fp_era20c_*.nc'],0,'WISC_era20c_eur_WS',20);
    else
        hazard_eraint=wisc_hazard_set([wisc_dir filesep 'fp_eraint_*.nc'],0,'WISC_eraint_eur_WS',20);
    end
    
end % loop_i
delete(pool)

hazard=climada_hazard_merge(hazard_era20c,hazard_eraint,'events');
hazard.filename = [hazards_dir filesep 'WISC_eur_WS' '.mat'];
fprintf('> saving combined hazard as %s\n',hazard.filename);
save(hazard.filename,'hazard',climada_global.save_file_version);

% copy results to dkrz (no closing ; to log success):
% ---------------------
%[status,result]=system('scp -r    /cluster/scratch/dbresch/climada_data/hazards/*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.')
%[status,result] =system('scp -r -v /cluster/scratch/dbresch/climada_data/hazards/*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.')

exit % the cluster appreciates this, gives back memory etc.