%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WISC Open Street Map Entity %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Thomas R??sli, thomas.roeoesli@usys.ethz.ch, 20171122, init
tic
% define location of WISC Tier 3 files (downloaded from WISC homepage)
wisc_dir = 'D:\Documents_DATA\WISC_data_20170918\Tier 3 Indicators';
exposure_dir = 'D:\Documents_DATA\WISC_Exposure\Exposure';
if ~isdir(wisc_dir)|~isdir(exposure_dir)
    % ask for folder
    sprints('Define wisc_dir and exposure_dir to show location of WISC_files.\n');
end
% country_names = {'Austria','Belgium','Czech Republic','Denmark', ...
%     'Finland','France','Germany','Ireland','Italy', ... 
%     'Luxembourg','Netherlands','Norway','Poland','Portugal', ...
%     'Spain','Sweden','United Kingdom','Estonia','Lithuania', ...
%     'Latvia','Switzerland'}; 
country_names = {'Austria','Belgium','Czech Republic','Denmark', ...
    'Finland','France','Germany','Ireland','Italy', ... 
    'Luxembourg','Netherlands','Poland','Portugal', ...
    'Spain','Sweden','United Kingdom','Estonia','Lithuania', ...
    'Latvia','Switzerland'}; % GDP share for Norway not available, hence deleted, GDP factors for Switzerland assumed 1
% country_names = {'Austria','Belgium','Czech Republic','Denmark', ...
%     'Finland','Germany','Ireland','Italy', ... 
%     'Luxembourg','Netherlands','Poland','Portugal', ...
%     'Spain','Sweden','United Kingdom','Estonia','Lithuania', ...
%     'Latvia'}; % without France, Norway, Switzerland (original files not good)
%load WISC entity (create it using wisc_demo.m first)
entity  = load([climada_global.entities_dir filesep 'entity_blackmarble_Europe_WISC_grid_20171113_encoded_hist.mat']);
entity = entity.entity_encoded_era20c;
entity.assets.Value(:) = 0;

% read reconstuction cost from file [in $ or euro per m^2]
[maximum_losses_temp, maximum_losses_ind] = xlsread([wisc_dir filesep 'C3S_WISC_Risk_and_Loss_Assessment_Inputs.xlsx'],'2. Maximum Losses per Country','A2:D26');

%read GDP per NUTS3
[NUTS_gdp_share, NUTS_name_gdp] = xlsread([wisc_dir filesep 'C3S_WISC_Risk_and_Loss_Assessment_Inputs.xlsx'],'3. GDP Ratio per NUTS3','A2:B1316');

% read CLC_Class
CLC_table = xlsread([wisc_dir filesep 'C3S_WISC_Risk_and_Loss_Assessment_Inputs.xlsx'],'4. Building Type per country_v1','K42:M85');
%CLC_table =[CLC_table; [0 0 0]]; % add zero reconstruction value if CLC_Class is 0
%% for country_i = DE
for country_i = 1:numel(country_names)
    [country_name_i, country_ISO3_i] = climada_country_name(country_names(country_i));
    country_ISO2_i = country_ISO3_i(1:2);
    if strcmp(country_ISO2_i,'AU'), country_ISO2_i = 'AT'; end
    if strcmp(country_ISO3_i,'GBR'), country_ISO2_i = 'UK'; end
    if strcmp(country_ISO3_i,'SWE'), country_ISO2_i = 'SE'; end
    if strcmp(country_ISO3_i,'DNK'), country_ISO2_i = 'DK'; end
    if strcmp(country_ISO3_i,'IRL'), country_ISO2_i = 'IE'; end
    if strcmp(country_ISO3_i,'POL'), country_ISO2_i = 'PL'; end
    if strcmp(country_ISO3_i,'PRT'), country_ISO2_i = 'PT'; end
    if strcmp(country_ISO3_i,'EST'), country_ISO2_i = 'EE'; end
    % read reconstuction cost from file [in $ or euro per m^2]
    maximum_losses_ind_country_i = find(strcmp(maximum_losses_ind, country_name_i));
    maximum_losses_country_i.residential = maximum_losses_temp(maximum_losses_ind_country_i,1);
    maximum_losses_country_i.commercial = maximum_losses_temp(maximum_losses_ind_country_i,2);
    maximum_losses_country_i.industrial = maximum_losses_temp(maximum_losses_ind_country_i,3);

    % WISC coordintes of country_i
    entity_ind_country_i = find(strcmp(entity.assets.centroid_admin0_ISO3,country_ISO3_i));
    entity_coord_country_i = [entity.assets.lat(entity_ind_country_i)' entity.assets.lon(entity_ind_country_i)'];

    % unpack file
    % untar([wisc_dir filesep 'C3S_WISC_RISK_' country_ISO2_i '_0100.tgz'], [wisc_dir filesep 'RISK']);
    
    % identify all NUTS3 files
    files_country_i = dir([exposure_dir filesep country_ISO2_i filesep '*_exposure.csv']);
    for file_i = 1:numel(files_country_i)
        fprintf(['reading file: ' files_country_i(file_i).name '.\n']);
        drawnow;
        %reading file
        data_temp = importdata([files_country_i(file_i).folder filesep files_country_i(file_i).name]);
        data_temp.('index') = str2double(data_temp.textdata(2:end,1));
        for i = 3:5
            data_temp.(data_temp.textdata{1,i}) = data_temp.textdata(2:end,i);
        end
        for i = 1:4
            data_temp.(data_temp.textdata{1,i+5}) = data_temp.data(:,i);
        end
        data_temp.CLC_2012(isnan(data_temp.CLC_2012)) = 0;
        data_temp = rmfield(data_temp,'data');
        data_temp = rmfield(data_temp,'textdata');

        %GDP factor on reconstruction costs
        GDP_factor_temp = NUTS_gdp_share(strcmp(NUTS_name_gdp,files_country_i(file_i).name(1:5)));
        if strcmp(country_ISO2_i, 'CH') 
            GDP_factor_temp = 1;
        end
%         if strcmp(country_ISO2_i, 'NO')
%             GDP_factor_temp = 1;
%         end
        % find centroid in entity closest to each entry
        next_neigbour_temp = knnsearch(entity_coord_country_i, [data_temp.LAT data_temp.LON]);

        % fill data into entity
        for building_i = 1:numel(data_temp.LAT)
            reconstruction_cost_i = CLC_table(CLC_table(:,1)==data_temp.CLC_2012(building_i),2) * ...
                mean([maximum_losses_country_i.residential maximum_losses_country_i.commercial]) ...
                + ...
                CLC_table(CLC_table(:,1)==data_temp.CLC_2012(building_i),3) * ...
                maximum_losses_country_i.industrial;
            if isempty(reconstruction_cost_i), reconstruction_cost_i = 0; end
            value_i = data_temp.AREA_m2(building_i) * reconstruction_cost_i * GDP_factor_temp;
            entity.assets.Value(entity_ind_country_i(next_neigbour_temp(building_i))) = ...
                value_i + entity.assets.Value(entity_ind_country_i(next_neigbour_temp(building_i)));
        end % building_i
    end % file_i
end % country_i
figure; climada_entity_plot(entity);
entity.assets.comment = 'Aggregated WISC reconstruction cost per Open Street Map Building';
entity.assets.filename = [climada_global.entities_dir filesep 'entity_WISC_OSM_Europe_WISC_grid_20171128_encoded_hist.mat'];
save(entity.assets.filename,'entity','-v7.3');
toc
% entity = climada_entity_load([climada_global.entities_dir filesep 'entity_WISC_OSM_Europe_WISC_grid_20171122_encoded_hist.mat']);
%% comparison
% entity_DEU = climada_entity_load([climada_global.entities_dir filesep 'DEU_Germany_blackmarble.mat']);
entity_Night = load([climada_global.entities_dir filesep 'entity_blackmarble_Europe_WISC_grid_20171113_encoded_hist.mat']);
entity_Night = entity_Night.entity_encoded_era20c;
figure; plot(entity_Night.assets.Value,entity.assets.Value, '.k')
% figure; climada_entity_plot(entity_DEU)
