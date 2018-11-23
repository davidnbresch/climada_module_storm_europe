%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% update of WISC analysis for revised synthetic datasets %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assessments of the revised WISC synthetic eventsets v1.2 v2 and v3 after the
% initally published synthetic eventset (v1 or original) showed too small
% storm systems (smaller areas affected compared to historic events) as
% shown in https://wisc.climate.copernicus.eu/wisc/documents/shared/(C3S_441_Lot3_WISC_SC2-D5.3-CGI-RP-18-0103)%20(ETH%20%20Swiss%20Re%20Case%20Study)%20(1.0).pdf


%% get the data of WISC
% get the synthetic eventset (1) as well as the historic footprints(2)
% (1) Please download the following files from the website:
% https://wisc.climate.copernicus.eu/wisc/#/help/products#eventset_section
% 
% 
% Summary files
% Containing maximum gust, mean gust and storm severity index (SSI) for each synthetic event:
% 
%     eventset_v1.2_summary.csv
%     eventset_v2_summary.csv
%     eventset_v3_summary.csv
% 
% Original WISC synthetic Event Set
% 
%     calibrated_v1_0001_to_1600.zip
%     calibrated_v1_1601_to_3200.zip
%     calibrated_v1_3201_to_4800.zip
%     calibrated_v1_4801_to_6400.zip
%     calibrated_v1_6401_to_7660.zip
% 
% WISC Synthetic Event Set 2
% 
%     calibrated_v2_0001_to_1600.zip
%     calibrated_v2_1601_to_3200.zip
%     calibrated_v2_3201_to_4800.zip
%     calibrated_v2_4801_to_6400.zip
%     calibrated_v2_6401_to_7660.zip
% 
% WISC Synthetic Event Set 3
% 
%     calibrated_v3_0001_to_1600.zip
%     calibrated_v3_1601_to_3200.zip
%     calibrated_v3_3201_to_4800.zip
%     calibrated_v3_4801_to_6400.zip
%     calibrated_v3_6401_to_7660.zip 
%
% further information here on the revised WISC event set:
% https://wisc.climate.copernicus.eu/wisc/documents/shared/(C3S_441_Lot3_WISC_SC2-D3.2.1-CGI-RP-17-0080)%20(C3S%20WISC%20Event%20Set%20)%20(1.1).pdf
synth_folder2 = 'D:\Documents_DATA\WISC_data_20180411\Event Set\'; % This folder should contain all the files mentioned above
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Please download the following files from the website:
% https://wisc.climate.copernicus.eu/wisc/#/help/products#footprint_section
% filename C3S_WISC_FOOTPRINT_NETCDF_0100.tgz
hist_folder = 'D:\Documents_DATA\WISC_data_20181029\Historic Storm Footprints\'; % This folder should contain the file mentioned above

%% quick look comparing the SSI


% synth_folder = 'D:\Documents_DATA\WISC_data_20170918\Synthetic Event Set\';
% SSI_synth_v1 = xlsread([synth_folder 'eventset_summary.csv'],'D2:D7661');


SSI_synth_v1_2 = xlsread([synth_folder2 'eventset_v1.2_summary.csv'],'D2:D7661');
SSI_synth_v2 = xlsread([synth_folder2 'eventset_v2_summary.csv'],'E2:E7661');
SSI_synth_v3 = xlsread([synth_folder2 'eventset_v3_summary.csv'],'E2:E7661');

% figure; boxplot([SSI_synth_v1 SSI_synth_v1_2 SSI_synth_v2 SSI_synth_v3],'labels',{'v1' 'v1.2' 'v2' 'v3'})
figure; boxplot([SSI_synth_v1_2 SSI_synth_v2 SSI_synth_v3],'labels',{'v1.2' 'v2' 'v3'})
ylabel('Storm Severity Index SSI')

% => Conclusion: SSI does indeed increase, which could mean that the area is indeed increasing 

%% quick look now comparing the area and mean gust speed of v2 and v3 to see if the increase in SSI in v3 can be attributed to bigger affected area
% using the reported area and mean gust speed of the WISC event set summary
% files
area_synth_v2 = xlsread([synth_folder2 'eventset_v2_summary.csv'],'B2:B7661');
area_synth_v3 = xlsread([synth_folder2 'eventset_v3_summary.csv'],'B2:B7661');

figure; boxplot([area_synth_v2 SSI_synth_v3],'labels',{'v2' 'v3'})
ylabel('area with gust speed > 25 m/s [km^2]')

gust_synth_v2 = xlsread([synth_folder2 'eventset_v2_summary.csv'],'D2:D7661');
gust_synth_v3 = xlsread([synth_folder2 'eventset_v3_summary.csv'],'D2:D7661');

figure; boxplot([gust_synth_v2 gust_synth_v3],'labels',{'v2' 'v3'})
ylabel('mean gust speed of gust speeds > 25 m/s [m/s]')

% => Conclusion: v3 has a much bigger affected area per storm event
% compared to v2, which could mean that the synthetic storms in v3 could be
% resembling the historic events much better v1.2 and v2


%% in depth look at area and mean gust speed of synthetic event set v1.2, v2 and v3 compared to the historic events

%create climada hazards based on the historic and synthetic events

% create historic hazard
hazard_wisc_hist_file=[climada_global.hazards_dir filesep 'WISC_eur_WS.mat'];
if exist('hazard_wisc_hist','var')
    fprintf('NOTE: using previously loaded hazard_wisc_hist\n');
elseif exist(hazard_wisc_hist_file,'file')
    fprintf('loading %s\n',hazard_wisc_hist_file);
    hazard_wisc_hist=climada_hazard_load(hazard_wisc_hist_file);
else % create WISC historic hazard
    wisc_hist_netcdf_glob = [hist_folder 'C3S_WISC_FOOTPRINT_NETCDF_0100' filesep 'fp_*.nc'];
    wisc_hist_netcdf_files = dir(wisc_hist_netcdf_glob);
    if numel(wisc_hist_netcdf_files) == 0
        % unzip downloaded zip file
        untar([hist_folder 'C3S_WISC_FOOTPRINT_NETCDF_0100.tgz'],[hist_folder 'C3S_WISC_FOOTPRINT_NETCDF_0100']);
        wisc_hist_netcdf_files = dir(wisc_hist_netcdf_glob);
        if numel(wisc_hist_netcdf_files) == 0, error(['ERROR: Historic footprints could not be '...
            'extracted from downloaded file:"C3S_WISC_FOOTPRINT_NETCDF_0100.tgz". '...
            'Please download this file or rename the variable "hist_folder" in '...
            'this script to match the file location or grant writing access '...
            'to this file location.']); 
        end
    end
    % generate the hazard event set from netCDF footprints
    % ====================================================
    hazard_wisc_hist=wisc_hazard_set(wisc_hist_netcdf_glob,0,hazard_wisc_hist_file);
end

% create synthetic hazard
wisc_synth_versions = {'v1' 'v2' 'v3'};
for version_i = 1:numel(wisc_synth_versions)
    hazard_wisc_synth_file=[climada_global.hazards_dir filesep 'WISC_eur_WS_synth_' wisc_synth_versions{version_i} '_mem5.mat'];
    if exist(hazard_wisc_synth_file,'file')
        fprintf(['all hazards for WISC synthetic event set version: ' wisc_synth_versions{version_i} 'seem to exist.\n']);
        %hazard_wisc_hist=climada_hazard_load(hazard_wisc_synth_file); % do not load, not necessary here
    else % create WISC synthetic hazard
        wisc_synth_netcdf_glob = [synth_folder2 'calibrated_' wisc_synth_versions{version_i} filesep 'fp_*.nc'];
        wisc_synth_netcdf_files = dir(wisc_synth_netcdf_glob);
        if numel(wisc_synth_netcdf_files) == 0
            % unzip downloaded zip file
            unzip([synth_folder2 'calibrated_' wisc_synth_versions{version_i} '_0001_to_1600.zip'],[synth_folder2 'calibrated_' wisc_synth_versions{version_i}]);
            unzip([synth_folder2 'calibrated_' wisc_synth_versions{version_i} '_1601_to_3200.zip'],[synth_folder2 'calibrated_' wisc_synth_versions{version_i}]);
            unzip([synth_folder2 'calibrated_' wisc_synth_versions{version_i} '_3201_to_4800.zip'],[synth_folder2 'calibrated_' wisc_synth_versions{version_i}]);
            unzip([synth_folder2 'calibrated_' wisc_synth_versions{version_i} '_4801_to_6400.zip'],[synth_folder2 'calibrated_' wisc_synth_versions{version_i}]);
            unzip([synth_folder2 'calibrated_' wisc_synth_versions{version_i} '_6401_to_7660.zip'],[synth_folder2 'calibrated_' wisc_synth_versions{version_i}]);
            wisc_synth_netcdf_files = dir(wisc_synth_netcdf_glob);
            if numel(wisc_synth_netcdf_files) == 0, error(['ERROR: synthetic footprints could not be '...
                'extracted from downloaded zip files. '...
                'Please download the files or rename the variable "synth_folder2" in '...
                'this script to match the file location or grant writing access '...
                'to this file location.']); 
            end
        end
        % generate the hazard event set from netCDF footprints
        % ====================================================
        for member_i = 1:5
            hazard_wisc_synth_mem_file = [climada_global.hazards_dir filesep 'WISC_eur_WS_synth_' ...
                wisc_synth_versions{version_i} '_mem' num2str(member_i) '.mat'];
            if exist(hazard_wisc_synth_mem_file,'file'), continue; end
            wisc_synth_mem_netcdf_glob = [synth_folder2 'calibrated_' ...
                wisc_synth_versions{version_i} filesep 'fp_*' num2str(member_i) '.nc'];
%             hazard_wisc_synth=wisc_hazard_set(wisc_synth_mem_netcdf_glob,0,hazard_wisc_synth_mem_file);
            hazard_wisc_synth=wisc_event_set_hazard(wisc_synth_mem_netcdf_glob,0,hazard_wisc_synth_mem_file);
            % maybe try wisc_event_set_hazard.m
            clear hazard_wisc_synth
        end
    end
end % version_i

% calculate the SSI, the affected area and the mean gust speed of the
% affected area for all event sets (hist and synthetic v1.2 v2 and v3)
% try using davids area calculation
if ~isfield(hazard_wisc_hist,'on_land')
    country_shapes=climada_shaperead(climada_global.map_border_file,1,1); % reads .mat
    hazard_wisc_hist.on_land = climada_inshape(hazard_wisc_hist.lat, hazard_wisc_hist.lon, country_shapes);
end
tic
[SSI_hist_th, SSI_struct_hist_th] = climada_hazard_SSI_calc(hazard_wisc_hist,25,hazard_wisc_hist.area_km2);
toc
tic
[SSI_hist_db, SSI_struct_hist_db] = climada_hazard_ssi(hazard_wisc_hist,1,25);
toc
% figure; plot(SSI_hist_th,SSI_hist_db,'xk') % this plot shows that the
% function "climada_hazard_SSI_calc" provides the same result as
% "climada_hazard_ssi", it is used here because it also outputs affected area and
% mean gust per event, which is important for 

% calc ssi for synthetic events
SSI_synth_db = cell(3,5);
SSI_struct_synth_db = cell(3,5);
SSI_synth_th = cell(3,5);
SSI_struct_synth_th = cell(3,5);
for version_i = 1:3
    for member_i = 1:5
        hazard_wisc_synth_file = [climada_global.hazards_dir filesep 'WISC_eur_WS_synth_v' ...
            num2str(version_i) '_mem' num2str(member_i) '.mat'];
        hazard_wisc_synth = climada_hazard_load(hazard_wisc_synth_file);
        tic
        if (max(hazard_wisc_hist.lat - hazard_wisc_synth.lat) < 10^-3) & (max(hazard_wisc_hist.lon - hazard_wisc_synth.lon) < 10^-3)
            hazard_wisc_synth.area_km2 = hazard_wisc_hist.area_km2;
            hazard_wisc_synth.on_land = hazard_wisc_hist.on_land;
        else
            tmp.lon = reshape(hazard_wisc_synth.lon,hazard_wisc_synth.lonlat_size);
            tmp.lat = reshape(hazard_wisc_synth.lat,hazard_wisc_synth.lonlat_size);
            tmp.area_km2=abs(tmp.lon(2:end,2:end)-tmp.lon(1:end-1,1:end-1))...
                .*cos(tmp.lat(1:end-1,1:end-1)./180.*pi)*111.12 .* ...
                abs(tmp.lat(2:end,2:end)-tmp.lat(1:end-1,1:end-1))*111.12;
            % fill last row/column
            tmp.area_km2(end+1,:)=tmp.area_km2(end,:);
            tmp.area_km2(:,end+1)=tmp.area_km2(:,end);
            hazard_wisc_synth.area_km2 = tmp.area_km2(:);
            clear tmp
            country_shapes=climada_shaperead(climada_global.map_border_file,1,1); % reads .mat
            hazard_wisc_synth.on_land = climada_inshape(hazard_wisc_synth.lat, hazard_wisc_synth.lon, country_shapes);

        end
        [SSI_synth_db{version_i,member_i}, SSI_struct_synth_db{version_i,member_i}] = climada_hazard_ssi(hazard_wisc_synth,1,25);
        toc
        [SSI_synth_th{version_i,member_i}, SSI_struct_synth_th{version_i,member_i}] = climada_hazard_SSI_calc(hazard_wisc_synth,25,hazard_wisc_synth.area_km2);

    end
end
%save('20181115_Wisc_new_synth_case.mat','-v7.3')
% Method David
SSI_synth_db_sorted = cell(3,1);
synth_area_th_sorted = cell(3,1);
synth_gust3_th_sorted = cell(3,1);
for version_i = 1:3
    SSI_synth_db_sorted{version_i}=[];
    synth_area_th_sorted{version_i}=[];
    synth_gust3_th_sorted{version_i}=[];
    for member_i = 1:5
        SSI_synth_db_sorted{version_i} = [SSI_synth_db_sorted{version_i}, SSI_synth_db{version_i,member_i}];
        synth_area_th_sorted{version_i} = [synth_area_th_sorted{version_i}, SSI_struct_synth_th{version_i,member_i}.area];
        synth_gust3_th_sorted{version_i} = [synth_gust3_th_sorted{version_i}, SSI_struct_synth_th{version_i,member_i}.gustspeed3];
    end
    SSI_synth_db_sorted{version_i} = sort(SSI_synth_db_sorted{version_i});
    synth_area_th_sorted{version_i} = sort(synth_area_th_sorted{version_i});
    synth_gust3_th_sorted{version_i} = sort(synth_gust3_th_sorted{version_i});
end
%save('20181115_Wisc_new_synth_case2.mat','-v7.3')

% plot
figure; boxplot([SSI_hist_db, SSI_synth_db_sorted{1:3}],...
    [ones(size(SSI_hist_db)), ones(size(SSI_synth_db_sorted{1}))*2,...
    ones(size(SSI_synth_db_sorted{2}))*3, ones(size(SSI_synth_db_sorted{3}))*4],...
    'Labels',{'WISC historic event set','WISC synthetic event set v1.2','WISC synthetic event set v2','WISC synthetic event set v3'})
title('Europe');ylabel('SSI');set(gcf,'color','white'); 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.66, 0.66]);

% figure; boxplot([SSI_hist_db, SSI_synth_db_sorted{version_i}(((end-(167*3))+1):end)], [ones(size(SSI_hist_db)), ones([167*3,1])'*2])

figure; boxplot([SSI_struct_hist_th.area, synth_area_th_sorted{1:3}],...
    [ones(size(SSI_hist_db)), ones(size(synth_area_th_sorted{1}))*2,...
    ones(size(synth_area_th_sorted{2}))*3, ones(size(synth_area_th_sorted{3}))*4],...
    'Labels',{'WISC historic event set','WISC synthetic event set v1.2','WISC synthetic event set v2','WISC synthetic event set v3'})
title('Europe');ylabel('Affected Area [km^2]');set(gcf,'color','white'); 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.66, 0.66]);

figure; boxplot([SSI_struct_hist_th.gustspeed3, synth_gust3_th_sorted{1:3}],...
    [ones(size(SSI_hist_db)), ones(size(synth_gust3_th_sorted{1}))*2,...
    ones(size(synth_gust3_th_sorted{2}))*3, ones(size(synth_gust3_th_sorted{3}))*4],...
    'Labels',{'WISC historic event set','WISC synthetic event set v1.2','WISC synthetic event set v2','WISC synthetic event set v3'})
title('Europe');ylabel('gustspeed^3 [(m/s)^3]');set(gcf,'color','white'); 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.66, 0.66]);

%% comparison of damages for Europe and Netherlands
% for v3 of the synthetic event set, compare the calculated damages with
% the historic event set
% create lfc graphs
% load hazard
hazard_wisc_hist = climada_hazard_load([climada_global.hazards_dir filesep 'WISC_eur_WS.mat']);
% load/create entity
entity_filename = 'entity_blackmarble_Europe_WISC_grid_20171113_encoded_hist.mat';
if exist([climada_global.entities_dir filesep entity_filename],'file')
    load([climada_global.entities_dir filesep entity_filename]);
else
    country_names = {'Austria','Belgium','Czech Republic','Denmark', ...
        'Finland','France','Germany','Greece','Hungary','Ireland','Italy', ... 
        'Luxembourg','Netherlands','Norway','Poland','Portugal','Slovakia', ...
        'Slovenia','Spain','Sweden','United Kingdom','Estonia','Lithuania', ...
        'Latvia','Switzerland'};
    parameters.entity_filename = entity_filename;
    entity_blackmarble_Europe_WISC_grid = climada_centroids_generate_blackmarble_entity(hazard_wisc_hist,country_names,parameters);
    entity_encoded_era20c = climada_assets_encode(entity_blackmarble_Europe_WISC_grid,hazard_era20c);
    save([climada_global.entities_dir filesep entity_filename],'entity_encoded_era20c','-v7.3');
end
% calculate damages historic footprints
EDS(1) = climada_EDS_calc(entity_encoded_era20c,hazard_wisc_hist);
YDS(1) = climada_EDS2YDS(EDS(1),hazard_wisc_hist);
figure; DFC_EDS(1) = climada_EDS_DFC(EDS(1))
figure; DFC_YDS(1) = climada_EDS_DFC(YDS(1))
% calculate damages synthetic event set version 3
for i = 1:5
    hazard = climada_hazard_load([climada_global.hazards_dir filesep 'WISC_eur_WS_synth_v3' ...
            '_mem' num2str(i) '.mat']);
    EDS_i(i) = climada_EDS_calc(entity_encoded_era20c,hazard);
    YDS_i(i) = climada_EDS2YDS(EDS_i(i),hazard);
%     figure; DFC_EDS(2) = climada_EDS_DFC(EDS(2))
%     figure; DFC_YDS(2) = climada_EDS_DFC(YDS(2))
end
% combine YDS_i
YDS(3).damage = [YDS_i(1:5).damage]
YDS(3).yyyy = 1:numel([YDS_i(1:5).yyyy])
YDS(3).frequency = [YDS_i(1:5).frequency]/5
YDS(3).hazard = YDS(2).hazard
YDS(3).assets = YDS(2).assets
YDS(3).annotation_name = 'WISC synth eur WS'
figure; climada_EDS_DFC(YDS([1 3]))
% figure; climada_EDS_DFC(YDS([1 3 4]))
% combine EDS_i
EDS(3).damage = [EDS_i(1:5).damage]
EDS(3).frequency = [EDS_i(1:5).frequency]/5
EDS(3).annotation_name = 'WISC synth eur WS'   
figure; climada_EDS_DFC(EDS([1 3]))


% load([climada_global.entities_dir filesep 'entity_blackmarble_Europe_WISC_grid_20171113_encoded_hist.mat']);
% hazard = climada_hazard_load([climada_global.hazards_dir filesep 'WS_Europe.mat']);
% hazard.fraction = zeros(size(hazard.intensity));
% hazard.fraction(hazard.intensity>0) = 1;
% entity_encoded_WS = climada_assets_encode(entity_encoded_era20c,hazard);
% EDS(4) = climada_EDS_calc(entity_encoded_WS,hazard);
% YDS(4) = climada_EDS2YDS(EDS(4),hazard);
% figure; DFC_EDS(4) = climada_EDS_DFC(EDS(4))
% figure; DFC_YDS(4) = climada_EDS_DFC(YDS(4))
% figure; climada_EDS_DFC(EDS([1 3 4]))
% figure; climada_EDS_DFC(YDS([1 3 4]))

figure; boxplot([EDS([1 3]).damage],repelem([0 1],[numel(EDS(1).damage),numel(EDS(3).damage)]),...
    'Labels',{'WISC historic event set','WISC synthetic event set'})
title('Europe'); ylabel('Event Damage [USD]'); set(gcf,'color','white'); 
figure; boxplot([YDS([1 3]).damage],repelem([0 1],[numel(YDS(1).damage),numel(YDS(3).damage)]),...
    'Labels',{'WISC historic event set','WISC synthetic event set'})
title('Europe'); ylabel('Annual Damage [USD]'); set(gcf,'color','white'); 
% figure; boxplot([YDS([1 3 4]).damage],repelem([0 1 2],[numel(YDS(1).damage),numel(YDS(3).damage),numel(YDS(4).damage)]),...
%     'Labels',{'WISC historic event set','WISC synthetic event set','climada event set'})
% title('Europe'); ylabel('Annual Damage [USD]'); set(gcf,'color','white'); 
save([climada_global.results_dir filesep '20181120_Workspace_WISC_Europe.mat'], 'EDS', 'EDS_i', 'YDS', 'YDS_i')
% load([climada_global.results_dir filesep 'Workspace_WISC_Europe.mat'])

% netherlands annual damages

entity_NLD_filename = 'NLD_Netherlands_blackmarble.mat';
if exist([climada_global.entities_dir filesep entity_NLD_filename],'file')
    entity = climada_entity_load(entity_NLD_filename); % select NLD_blackmarble....mat
else
    parameters.entity_filename = entity_NLD_filename;
    entity = climada_centroids_generate_blackmarble_entity(hazard_wisc_hist,{'Netherlands'},parameters);
end

EDS(1) = climada_EDS_calc(entity,hazard_wisc_hist);
YDS(1) = climada_EDS2YDS(EDS(1),hazard_wisc_hist);
figure; DFC_EDS(1) = climada_EDS_DFC(EDS(1))
figure; DFC_YDS(1) = climada_EDS_DFC(YDS(1))

for i = 1:5
    hazard = climada_hazard_load([climada_global.hazards_dir filesep 'WISC_eur_WS_synth_v3' ...
        '_mem' num2str(i) '.mat']);
    if i==1, entity = climada_assets_encode(entity,hazard); end
    EDS_i(i) = climada_EDS_calc(entity,hazard);
    YDS_i(i) = climada_EDS2YDS(EDS_i(i),hazard);
%     figure; DFC_EDS(2) = climada_EDS_DFC(EDS(2))
%     figure; DFC_YDS(2) = climada_EDS_DFC(YDS(2))
end
% combine YDS_i
YDS(3).damage = [YDS_i(1:5).damage]
YDS(3).yyyy = 1:numel([YDS_i(1:5).yyyy])
YDS(3).frequency = [YDS_i(1:5).frequency]/5
YDS(3).hazard = YDS(2).hazard
YDS(3).assets = YDS(2).assets
YDS(3).annotation_name = 'WISC synth eur WS'
figure; climada_EDS_DFC(YDS([1 3]))
% combine EDS_i
EDS(3).damage = [EDS_i(1:5).damage]
EDS(3).frequency = [EDS_i(1:5).frequency]/5
EDS(3).annotation_name = 'WISC synth eur WS'
%EDS(3).yyyy = [EDS_i(1:5).yyyy]-repelem([1985 (1985-26) (1985-2*26) (1985-3*26) (1985-4*26)],...
    
figure; climada_EDS_DFC(EDS([1 3]))

% hazard = climada_hazard_load([climada_global.hazards_dir filesep 'WS_Europe.mat']);
% hazard.fraction = sparse(zeros(size(hazard.intensity)));
% hazard.fraction(hazard.intensity>0) = 1;
% entity = climada_assets_encode(entity,hazard);
% EDS(4) = climada_EDS_calc(entity,hazard);
% YDS(4) = climada_EDS2YDS(EDS(4),hazard);
% figure; DFC_EDS(4) = climada_EDS_DFC(EDS(4))
% figure; DFC_YDS(4) = climada_EDS_DFC(YDS(4))

figure; boxplot([EDS([1 3]).damage],repelem([0 1],[numel(EDS(1).damage),numel(EDS(3).damage)]),...
    'Labels',{'WISC historic event set','WISC synthetic event set'})
title('Netherlands'); ylabel('Event Damage [USD]'); set(gcf,'color','white'); 
figure; boxplot([YDS([1 3]).damage],repelem([0 1],[numel(YDS(1).damage),numel(YDS(3).damage)]),...
    'Labels',{'WISC historic event set','WISC synthetic event set'})
title('Netherlands'); ylabel('Annual Damage [USD]'); set(gcf,'color','white'); 
% figure; boxplot([YDS([1 3 4]).damage],repelem([0 1 2],[numel(YDS(1).damage),numel(YDS(3).damage),numel(YDS(4).damage)]),...
%     'Labels',{'WISC historic event set','WISC synthetic event set','climada event set'})
% title('Netherlands'); ylabel('Annual Damage [USD]'); set(gcf,'color','white'); 
save([climada_global.results_dir filesep '20181122_Workspace_WISC_Netherlands.mat'], 'EDS', 'EDS_i', 'YDS', 'YDS_i')


