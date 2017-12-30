function [intensity_prob,n_prob_events]=climada_ws_hist2prob(intensity2d,lonlat_size)
% climada historic probabilistic
% MODULE:
%   storm_europe
% NAME:
%   climada_ws_hist2prob
% PURPOSE:
%   Given a (historic) single event wind footprint, generate the probabilistic
%   events. Generation of probabilistic (daughter) events similar as
%   in Schwierz et al, 2010. See also old catXos code catxos_WSEU_generate_probab_set.m
%   and possibly the very old process_allme_probab.m and me2probailistic.m
%   (available upon request, not part of the open source distribution)
%
%   previous call: generate a winter storm hazard event set, e.g. see wisc_hazard_set
%   next call: check using e.g. climada_hazard_ssi
% CALLING SEQUENCE:
%   [intensity_prob,n_prob_events]=climada_ws_hist2prob(intensity2d,lonlat_size)
% EXAMPLE:
%   hazard=wisc_hazard_set('test'); % create a single event windstorm hazard
%   intensity_prob=climada_ws_hist2prob(hazard.intensity,hazard.lonlat_size)
% INPUTS:
%   intensity2d: a single hazard footprint, either as 2d or 1d field. In
%       case of 1d, lonlat_size needs to be provided. We allow passing a 1d
%       array, as climada stores the intensity as (1,n_centroids) and hence
%       it is sometimes convenient to pass on that form. Passing on 2d is
%       faster, though.
% OPTIONAL INPUT PARAMETERS:
%   lonlat_size: the size of intensity2d, only required if intensity2d is
%       passed on as 1d array.
% OUTPUTS:
%   intensity_prob(n_prob_events+1,n_centroids): the probabilistic
%       intensity, first event is the original footprint, NOT sparse
%   n_prob_events: the number of probabilistic events added
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20171229, initial
% David N. Bresch, david.bresch@gmail.com, 20171230, NOT sparse returned
%-

intensity_prob=[];

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('intensity2d','var'),return;end
if ~exist('lonlat_size','var'),lonlat_size=[];end

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
n_prob_events=20;
dii=4; % shift North/South and East/West by dii gridpoints (use e.g. 100 for checks)

if isempty(lonlat_size),lonlat_size=size(intensity2d);end
n_centroids=prod(lonlat_size);
intensity2d0=zeros(lonlat_size); % init empty array
% faster (facor two) to first deal with double arrays, then sparsify
%intensity_prob=spalloc(1+n_prob_events,n_centroids,(1+n_prob_events)*n_centroids*nnz(intensity2d)/numel(intensity2d));
intensity_prob=zeros(1+n_prob_events,n_centroids);

if size(intensity2d,1)==1
    intensity_prob(1,:)=intensity2d;next_event=2;
    if ~isempty(lonlat_size)
        % convert back to 2d-array
        %lon2d=reshape(hazard.lon,hazard.lonlat_size);
        %lat2d=reshape(hazard.lat,hazard.lonlat_size);
        intensity2d=reshape(intensity2d(1,:),lonlat_size); % check with: contour(lon2d,lat2d,intensity2d)
    else
        fprintf('ERROR: either provide 2d intensity or lonlat_size in addition\n');
        return
    end
else
    intensity_prob(1,:)=reshape(intensity2d,1,n_centroids);next_event=2;
end

intensity2d_sqrt = sqrt(intensity2d); % for speed-up
intensity2d_pwrd =      intensity2d.^1.1; % for speed-up

for loop_i=1:5
    
    switch loop_i
        case 1
            intensity2d_modi=intensity2d; % translation only
        case 2
            %intensity2d_modi=intensity2d-0.06*intensity2d_sqrt;
            intensity2d_modi=intensity2d-0.1*intensity2d_sqrt;
        case 3
            %intensity2d_modi=intensity2d+0.15*intensity2d_sqrt;
            intensity2d_modi=intensity2d+0.1*intensity2d_sqrt;
        case 4
            %         intensity2d_modi=intensity2d_temp-0.0040*intensity2d_temp;
            %         intensity2d_modi=intensity2d_temp+0.0125*intensity2d_temp;
            %         intensity2d_modi=intensity2d_temp+0.1200*intensity2d_temp;
            %intensity2d_modi=intensity2d_temp+0.04*intensity2d_pwrd;
            intensity2d_modi=intensity2d+0.1*intensity2d_pwrd;
        case 5
            %intensity2d_modi=intensity2d_temp-0.0330*intensity2d_temp;
            %intensity2d_modi=intensity2d_temp-0.05*intensity2d_pwrd;
            intensity2d_modi=intensity2d-0.1*intensity2d_pwrd;
            %intensity2d_modi=intensity2d_temp-0.1200*intensity2d_temp;
    end % switch
    
    % shift North/South and East/West dii gridpoints
    intensity2d_temp=intensity2d0; % init
    intensity2d_temp(:,1:end-dii)=intensity2d_modi(:,1+dii:end); % South
    intensity_prob(next_event,:)=reshape(intensity2d_temp,1,n_centroids);next_event=next_event+1;
    intensity2d_temp=intensity2d0; % init
    intensity2d_temp(:,1+dii:end)=intensity2d_modi(:,1:end-dii); % North
    intensity_prob(next_event,:)=reshape(intensity2d_temp,1,n_centroids);next_event=next_event+1;
    
    intensity2d_temp=intensity2d0; % init
    intensity2d_temp(1:end-dii,:)=intensity2d_modi(1+dii:end,:); % West
    intensity_prob(next_event,:)=reshape(intensity2d_temp,1,n_centroids);next_event=next_event+1;
    intensity2d_temp=intensity2d0; % init
    intensity2d_temp(1+dii:end,:)=intensity2d_modi(1:end-dii,:); % East
    intensity_prob(next_event,:)=reshape(intensity2d_temp,1,n_centroids);next_event=next_event+1;
    
end % loop_i

%intensity_prob=sparse(intensity_prob); % sparsify

end % climada_ws_hist2prob