%function ZAMG_test
% climada test for ZAMG
% MODULE:
%   storm europe
% NAME:
%   ZAMG_test
% PURPOSE:
%   batch job to TEST ZAMG damage calculations
%
% CALLING SEQUENCE:
%   ZAMG_test
% EXAMPLE:
%   ZAMG_test
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to stdout and figures and to a test testdata.csv
%   results folder
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20180702, initial
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
%
% set year, month and day for the single event to treat
yyyy=1999;mm=12;dd=26;
%
hazard_name='WISC_AUT_eur_WS';
entity_name='AUT_Austria_01x01'; % needs to start with ISO3 country code
%entity_name='AUT_Austria_10x10';
%
markersize=2; % size of dots in plots
%
fig_dir = [climada_global.results_dir filesep 'ZAMG'];
if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
fig_ext ='png';

% prepare the asset base
entity_file=[climada_global.entities_dir filesep entity_name '.mat'];
if exist(entity_file,'file')
    entity=climada_entity_load(entity_file); % loads asset base
else
    entity=climada_entity_counbtry(entity_name(1:3),1); % create asset base
end

% prepare the hazard
hazard_file=[climada_global.hazards_dir filesep hazard_name  '.mat'];
if exist(hazard_file,'file')
    hazard=climada_hazard_load(hazard_file); % loads Copernicus WISC hazard event set for Austria
else
    hazard=wisc_hazard_set;
end

% figure position of the single event
event_index=find(hazard.yyyy==yyyy & hazard.mm==mm & hazard.dd==dd & hazard.orig_event_flag==1);
if length(event_index)==1
else
    fprintf('ERROR: no single storm found for %i%i%i, aborted\n',yyyy,mm,dd);
    return
end

% matches assets with hazard resolution
entity=climada_assets_encode(entity,hazard); 

% create a dummy hazard with only one single event
hazard_1=hazard; % copy, now only keep one event:
hazard_1.intensity       = hazard_1.intensity(event_index,:);
hazard_1.fraction        = hazard_1.fraction(event_index,:);
hazard_1.frequency       = hazard_1.frequency(event_index);
hazard_1.event_ID        = hazard_1.event_ID(event_index);hazard_1.yyyy=hazard_1.yyyy(event_index);
hazard_1.orig_event_flag = hazard_1.orig_event_flag(event_index);
hazard_1.mm=hazard_1.mm(event_index);hazard_1.dd=hazard_1.dd(event_index);

% calculate damage for single event
EDS=climada_EDS_calc(entity,hazard_1,'single');

fprintf('%s %2.2f million\n',EDS.Value_unit,EDS.damage/1e6); % the uncalibrated Lothar damage in USD)

% and now the plots
figure;climada_entity_plot(entity,markersize); % plot the asset base
title('Austria - exposure (1x1km resolution)')
saveas(gcf,[fig_dir filesep 'ZAMG_assets'],fig_ext);

figure;climada_hazard_plot(hazard,event_index,markersize); % plots Lothar (happens to be at position event_index in the fully probabilistic hazard set)
title('Austria - Lothar wind field (4x4km resolution)')
saveas(gcf,[fig_dir filesep 'ZAMG_hazard'],fig_ext);

damage_data=entity;damage_data.assets.Value=EDS.ED_at_centroid'/hazard_1.frequency;
fprintf('check: %s %2.2f million\n',EDS.Value_unit,sum(damage_data.assets.Value)/1e6); % the uncalibrated Lothar damage in USD)
figure;climada_entity_plot(damage_data,markersize);
title('Austria - Simulated Lothar damage (1x1km resolution)')
saveas(gcf,[fig_dir filesep 'ZAMG_damage'],fig_ext);

% just store in one variable for ease of readinf write statement below
data.lon          = entity.assets.lon; 
data.lat          = entity.assets.lat;
data.asset_value  = entity.assets.Value;
data.damage_value = damage_data.assets.Value;

fprintf('exporting testdata ...');
fid=fopen('testdata.csv','w'); % open raw text file
for i=1:length(data.lon)
   fprintf(fid,'%f;%f;%f;%f\n',data.lon(i),data.lat(i),data.asset_value(i),data.damage_value(i));
end
fclose(fid);
fprintf(' done\n');
