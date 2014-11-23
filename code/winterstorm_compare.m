function EDS=winterstorm_compare(entity)
% climada
% NAME:
%   winterstorm_compare
% PURPOSE:
%   Run a winterstorm Europe analysis for a given entity with different
%   damage functions (eg those derived based upon DWD and EuroTempest gust
%   data). The code tries to use the matching hazard event set (as as
%   created by country_risk_calc). If not found, the code prompts for a WS
%   hazard event set and re-encodes assets.
%
% CALLING SEQUENCE:
%   EDS=winterstorm_compare(entity)
% EXAMPLE:
%   EDS=winterstorm_compare
% INPUTS:
%   entity: an encoded entity, see climada_entity_read
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   EDS: the event damage set(s), see climada_EDS_calc
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141121, ICE initial
%-

EDS=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('entity','var'),entity=[];end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
%
unique_DamageFunID=1234; % a kind of unique ID, unlikely to exist
%
%
storm_info.name={'Daria' 'Lothar' 'Kyrill'};

% Storm	Date*	Insured loss (USD, indexed to 2012)+	Affected countries	Umax**(ms-1)	Lowest MSLP++ (hPa)	Maximum vorticity***(10-5-1)	Sft+++
% Daria (Burns' Day storm)	25/1/1990	8.2bn	Belgium, France, Germany, Netherlands and United Kingdom	37.92	948.62	11.89	48.05
% Lothar	26/12/1999	8.0bn	France, Germany and Switzerland	36.72	973.16	6.43	18.82
% Kyrill	18/1/2007	6.7bn	Austria, Belgium, France, Germany, Ireland, Netherlands and United Kingdom	36.38	961.12	8.55	59.43
% Great Storm of 87	16/10/1987	6.3bn	France and United Kingdom	39.53	955.97	10.17	38.42
% Vivian	26/2/1990	5.6bn	Belgium, France, Germany, Netherlands and United Kingdom	35.16	940.30	9.53	40.86
% Klaus	24/1/2009	3.5bn	Andorra, France, Germany, Italy, Spain and Switzerland	37.23	965.95	9.25	24.36
% Martin	27/12/1999	3.3bn	France, Italy and Switzerland	37.18	967.58	9.18	21.33
% Xynthia	27/2/2010	2.9bn	Belgium, Denmark, France, Germany, Poland, Portugal, Spain, Sweden and United Kingdom	32.62	965.46	9.58	23.11
% Anatol	3/12/1999	2.6bn	Denmark, Germany and Sweden	39.86	956.05	10.98	47.01
% Erwin (Gudrun)	8/1/2005	2.2bn	Denmark, Ireland, Norway, Sweden and United Kingdom	39.22	959.89	9.82	36.08
% Herta	3/2/1990	1.5bn	Belgium, France, Germany, Netherlands and United Kingdom	33.16	944.78	13.36	15.94
% Emma	29/2/2008	1.4bn	Austria, Belgium, Czech Repubic, Germany, Netherlands, Poland and Switzerland	25.12	959.46	9.60	12.17
% Wiebke	28/2/1990	1.4bn	Belgium, France, Germany, Netherlands, Switzerland and United Kingdom	32.24	944.18	7.77	25.16
% Gero	11/1/2005	0.6bn	Ireland and United Kingdom	39.13	960.64	8.55	17.55
% Ulli	3/1/2012	0.2bn	United Kingdom	36.32	954.30	10.24	19.02
% Dagmar (Patrick)	26/12/2011	0.04bn	Finland and Norway	30.08	953.94	8.58	1.77    

% Storm	Date*	Insured loss (USD, indexed to 2012)+	Affected countries	Umax**(ms-1)	Lowest MSLP++ (hPa)	Maximum vorticity***(10-5-1)	Sft+++
% Daria (Burns' Day storm)	25/1/1990	8.2bn	Belgium, France, Germany, Netherlands and United Kingdom	37.92	948.62	11.89	48.05
% Lothar	26/12/1999	8.0bn	France, Germany and Switzerland	36.72	973.16	6.43	18.82
% Kyrill	18/1/2007	6.7bn	Austria, Belgium, France, Germany, Ireland, Netherlands and United Kingdom	36.38	961.12	8.55	59.43
% Great Storm of 87	16/10/1987	6.3bn	France and United Kingdom	39.53	955.97	10.17	38.42
% Vivian	26/2/1990	5.6bn	Belgium, France, Germany, Netherlands and United Kingdom	35.16	940.30	9.53	40.86
% Klaus	24/1/2009	3.5bn	Andorra, France, Germany, Italy, Spain and Switzerland	37.23	965.95	9.25	24.36
% Martin	27/12/1999	3.3bn	France, Italy and Switzerland	37.18	967.58	9.18	21.33
% Xynthia	27/2/2010	2.9bn	Belgium, Denmark, France, Germany, Poland, Portugal, Spain, Sweden and United Kingdom	32.62	965.46	9.58	23.11
% Anatol	3/12/1999	2.6bn	Denmark, Germany and Sweden	39.86	956.05	10.98	47.01
% Erwin (Gudrun)	8/1/2005	2.2bn	Denmark, Ireland, Norway, Sweden and United Kingdom	39.22	959.89	9.82	36.08
% Herta	3/2/1990	1.5bn	Belgium, France, Germany, Netherlands and United Kingdom	33.16	944.78	13.36	15.94
% Emma	29/2/2008	1.4bn	Austria, Belgium, Czech Repubic, Germany, Netherlands, Poland and Switzerland	25.12	959.46	9.60	12.17
% Wiebke	28/2/1990	1.4bn	Belgium, France, Germany, Netherlands, Switzerland and United Kingdom	32.24	944.18	7.77	25.16
% Gero	11/1/2005	0.6bn	Ireland and United Kingdom	39.13	960.64	8.55	17.55
% Ulli	3/1/2012	0.2bn	United Kingdom	36.32	954.30	10.24	19.02
% Dagmar (Patrick)	26/12/2011	0.04bn	Finland and Norway	30.08	953.94	8.58	1.77

% prompt for entity if not given
if isempty(entity),entity=climada_entity_load;end
if isempty(entity),return;end

% We use the WS Europe hazard event set as created by country_risk_calc
%
% Note that one would need to re-encode assets to each hazard prior to
% calling the damage calculation (as climada_EDS_calc assumes matching
% order for speedup), unless one knows that all hazard event sets are valid
% on the exact same centroids (the first n elements in the hazard are
% matching the n locations of the assets, while the n+1:end elements in
% hazard are the ones for the buffer around the country). The call to
% centroids_generate_hazard_sets ensures that, and hence no need for
% re-encoding (force_re_encoding=0). But in case you're in doubt, set
% force_re_encoding=1 and check whether you get the same results (if yes,
% very likely no need for force_re_encoding=1, otherwise keep =1).
% Note further that the code checks re-encoding (issuing a WARNING) if it
% looks like assets and hazard do not match and then sets
% force_re_encoding=1
force_re_encoding=0; % default=0

WS_hazard_event_set_file=[climada_global.data_dir filesep 'hazards' filesep deblank(entity.assets.filename) '_WS_Europe.mat'];

if ~exist(WS_hazard_event_set_file,'file')
    % propt for WS hazard event set
    WS_hazard_event_set_file=[module_data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(WS_hazard_event_set_file, 'Select a WS hazard event set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        WS_hazard_event_set_file=fullfile(pathname,filename);
        force_re_encoding=1; % to be on the safe side
    end
end

load(WS_hazard_event_set_file); % loads hazard

% check for entity to be encoded to the hazard event set's centroids
% in essence, re-encoce if in doubt (as otherwise results might be useless)
n_centroids=length(entity.assets.Longitude);
if size(hazard.lon,1)>1,hazard.lon=hazard.lon';hazard.lat=hazard.lat';end % transpose old hazard event sets
if length(hazard.lon)>=n_centroids
    if sum(abs(entity.assets.Longitude-hazard.lon(1:n_centroids)))+...
            sum(abs(entity.assets.Latitude-hazard.lat(1:n_centroids)))>0
        fprintf('WARNING: assets might not be encoded to hazard centroids\n');
        force_re_encoding=1;
    end
else
    fprintf('WARNING: assets not encoded to hazard centroids\n');
    force_re_encoding=1; % sure to re-encode, not even vector lengths match
end

if force_re_encoding
    fprintf('re-encoding hazard to the respective centroids\n');
    assets = climada_assets_encode(entity.assets,hazard);
    entity=rmfield(entity,'assets');
    entity.assets=assets; % assign re-encoded assets
end

EDS=climada_EDS_calc(entity,hazard);
EDS.annotation_name='Default';
next_EDS=2;

for test_i=1:2
    
    if test_i==1
        % with damage function as derived based on DWD gust data
        MDR_exp=0.2539;
        annotation_name='DWD';
    elseif test_i==2
        % with damage function as derived based on Euro Tempest gust data
        MDR_exp=0.2493;
        annotation_name='EuroTempest';
    end
    
    % calculate the damage function
    entity.assets.DamageFunID=entity.assets.DamageFunID*0+unique_DamageFunID;
    Intensity=0:2:100;
    MDR=0.00000001*exp(MDR_exp*Intensity);
    MDR=min(MDR,1);
    
    % switch damage function
    entity.damagefunctions.DamageFunID=Intensity*0+unique_DamageFunID;
    entity.damagefunctions.Intensity=Intensity;
    entity.damagefunctions.MDD=MDR;
    entity.damagefunctions.PAA=Intensity*0+1;
    entity.damagefunctions.MDR=MDR;
    if isfield(entity.damagefunctions,'peril_ID'),entity.damagefunctions=rmfield(entity.damagefunctions,'peril_ID');end
    
    % calculate EDS
    EDS(next_EDS)=climada_EDS_calc(entity,hazard);
    EDS(next_EDS).annotation_name=annotation_name;
    next_EDS=next_EDS+1;
    
end % test_i

climada_EDS_DFC(EDS,'',1);

return
