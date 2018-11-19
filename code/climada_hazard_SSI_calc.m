function [ SSI,SSI_struct ] = climada_hazard_SSI_calc(hazard,gustspeed_threshold,area_per_centroid)
%climada_hazard_SSI_calc - calculates the Dawkins (2016) Definition of
%Storm Severity Index (SSI) for a climada hazard containing winterstorm
%events
%   Detailed explanation goes here
global climada_global
if ~exist('area_per_centroid','var'),area_per_centroid=[];end
if ~exist('gustspeed_threshold','var')
    gustspeed_threshold=25;
    fprintf('gustspeed threshold not provided. 25 m/s used. unit of intensity in hazard is assumed to be m/s. \n');
end
if ~exist('hazard','var'), return; end

if isempty(area_per_centroid)
    fprintf('Area per centroid not provided: Estimate it based on 4 nearest centroids.\n');
    [~, min_dist] = knnsearch([hazard.lat;hazard.lon]',...
        [hazard.lat;hazard.lon]','k',4);
    km_per_degree = 110*cosd(hazard.lat); % estimating the distance per degree (VERY BAD implementation of distance)
    area_per_centroid = ((mean(min_dist,2)/2).*km_per_degree') .^2; 
end

if size(area_per_centroid) ~= size(hazard.intensity(1,:))
    area_per_centroid = reshape(area_per_centroid,size(hazard.intensity(1,:)));
end

% cut to land area
if isfield(hazard,'on_land')
    land_sea_mask = hazard.on_land;
else
    country_shapes=climada_shaperead(climada_global.map_border_file,1,1); % reads .mat
    land_sea_mask = climada_inshape(hazard.lat, hazard.lon, country_shapes);
    clear country_shapes;
end

% calculate SSI for each storm
SSI_struct.area = zeros(size(hazard.event_ID));
SSI_struct.gustspeed = zeros(size(hazard.event_ID));
SSI_struct.gustspeed3 = zeros(size(hazard.event_ID));
SSI_struct.SSI = zeros(size(hazard.event_ID));
for event_i = 1:hazard.event_count
    sel_pos = hazard.intensity(event_i,:)>gustspeed_threshold & ...
        land_sea_mask;
    
    SSI_struct.area(event_i) = sum(area_per_centroid(sel_pos));
    % SSI_struct.gustspeed(event_i) = mean(hazard.intensity(event_i,sel_pos));
    SSI_struct.gustspeed(event_i) = sum( ...
        hazard.intensity(event_i,sel_pos) .* area_per_centroid(sel_pos)) ...
        / SSI_struct.area(event_i);
    SSI_struct.gustspeed3(event_i) = SSI_struct.gustspeed(event_i)^3;
    SSI_struct.SSI(event_i) = SSI_struct.area(event_i) * SSI_struct.gustspeed3(event_i);
end

SSI_struct.area_per_centroid = area_per_centroid;
SSI_struct.land_sea_mask = land_sea_mask;

SSI = SSI_struct.SSI;
end

