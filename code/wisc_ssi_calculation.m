%% crosscheck of SSI and the involved mean gust speed and affected area between CLIMADA method and reported numbers by wisc
% in the storm europe module the WISC dataset uses the storm severity index SSI. 
% CLIMADA is calculating the SSI index as well (ssi = climada_hazard_ssi(hazard) )
% The reported values of the SSI differ from the CLIMADA reported values.
% Here we explain why and how.
% Please also read WISC_ssi_comparison.pdf in the folder "docs" of the  

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
synth_folder2 = 'D:\Documents_DATA\WISC_data_20180411\Event Set\'; % This folder should contain all the files mentioned above
% (2) Please create the hazard files from the WISC synthetic event set as
% done in the script WISC_new_synthetic_sets_case_study.m
if exist('hazard_wisc_synth','var') && str2double(hazard_wisc_synth.filename([end-9 end-4])) == 35
    % correct hazard still loaded
else
    fprintf('select hazard WISC_eur_WS_synth_v3_mem5.mat manually...\n')
    hazard_wisc_synth = climada_hazard_load();
    fprintf('select hazard WISC_eur_WS.mat manually...\n')
    hazard_wisc_hist = climada_hazard_load();
    hazard_wisc_synth.on_land = hazard_wisc_hist.on_land;
    hazard_wisc_synth.area_km2 = hazard_wisc_hist.area_km2;
end


% read SSI values as reported by WISC
SSI_synth_v3 = xlsread([synth_folder2 'eventset_v3_summary.csv'],'E2:E7661');
gust_synth_v3 = xlsread([synth_folder2 'eventset_v3_summary.csv'],'D2:D7661');
area_synth_v3 = xlsread([synth_folder2 'eventset_v3_summary.csv'],'B2:B7661');

version_i = 3; member_i = 5;
[SSI_synth_db{version_i,member_i}, SSI_struct_synth_db{version_i,member_i}] = climada_hazard_ssi(hazard_wisc_synth,0,25);
[SSI_synth_th{version_i,member_i}, SSI_struct_synth_th{version_i,member_i}] = climada_hazard_SSI_calc(hazard_wisc_synth,25,hazard_wisc_synth.area_km2);

% sort the reported numbers per member as is done in our process above
[~,filenames_synth_v3] = xlsread([synth_folder2 'eventset_v3_summary.csv'],'A2:A7661');
[~,order_filenames_synth_v3] = sort(cellfun(@(x) str2double(x(end-3)),filenames_synth_v3));
SSI_synth_v3_mem5 = SSI_synth_v3(order_filenames_synth_v3((end-hazard_wisc_synth.event_count+1):end));
area_synth_v3_mem5 = area_synth_v3(order_filenames_synth_v3((end-hazard_wisc_synth.event_count+1):end));
gust_synth_v3_mem5 = gust_synth_v3(order_filenames_synth_v3((end-hazard_wisc_synth.event_count+1):end));

% plot WISC reported mean gust vs. calculated mean gust
figure; plot(SSI_synth_v3_mem5,SSI_synth_th{3,5},'or');
xlabel('SSI, reported by WISC '); ylabel('SSI, calculated by CLIMADA '); 

% plot WISC reported area vs. calculated area
figure; plot(area_synth_v3_mem5,SSI_struct_synth_th{3,5}.area,'ob')
xlabel('Affected Area, reported by WISC [km^2]'); ylabel('Affected Area, calculated by CLIMADA [km^2]'); 

% plot WISC reported mean gust vs. calculated mean gust
figure; plot(gust_synth_v3_mem5,SSI_struct_synth_th{3,5}.gustspeed,'og')
xlabel('Mean gustspeed, reported by WISC [m/s]'); ylabel('Affected Area, calculated by CLIMADA [m/s]'); 

% area and mean gust are both quite different.... why?
% is SSI calulated by multiplying area with gust^3
all(abs(SSI_synth_v3-(area_synth_v3.*(gust_synth_v3).^3))<3) % true
% => yes SSI is calculated that way. 
% But then why are area and mean gust speeds so different....
% hypothesis one: WISC uses unit area
% test: use unit area as well:
[SSI_synth_th_test_unitarea, SSI_struct_synth_th_test_unitarea] = climada_hazard_SSI_calc(hazard_wisc_synth,25,ones(size(hazard_wisc_synth.area_km2))*11.69);
figure; plot(area_synth_v3_mem5,SSI_struct_synth_th_test_unitarea.area,'xk')
xlabel('Affected Area, reported by WISC [km^2]'); ylabel('Affected Area, calculated by CLIMADA [km^2]'); 
title('hypothesis 1: WISC uses unit area')
% => unit area does not explain the difference in area, 

figure; plot(SSI_synth_v3_mem5,SSI_synth_th{3,5},'oy'); hold on
plot(SSI_synth_v3_mem5,SSI_synth_th_test_unitarea,'xk')
% => SSI also still unmatching, our is always equal or bigger, maybe
% because we include iceland and finland and russia?
% hypothesis two: WISC uses different land_sea_mask
% test: use landmask without iceland and russia as well
hazard_wisc_synth_no_russia_iceland = hazard_wisc_synth; 
[country_name,~,~,incountry_russia] = climada_country_name('Russia','',hazard_wisc_synth_no_russia_iceland.lon,hazard_wisc_synth_no_russia_iceland.lat);
[country_name,~,~,incountry_iceland] = climada_country_name('Iceland','',hazard_wisc_synth_no_russia_iceland.lon,hazard_wisc_synth_no_russia_iceland.lat);
hazard_wisc_synth_no_russia_iceland.on_land(incountry_russia | incountry_iceland) = false;

[SSI_synth_th_noRUSICE, SSI_struct_synth_th_noRUSICE] = climada_hazard_SSI_calc(hazard_wisc_synth_no_russia_iceland,25,hazard_wisc_synth_no_russia_iceland.area_km2);
% plot old (color) and new (black) SSI, area, gust vs WISC reported numbers
figure; plot(area_synth_v3_mem5, SSI_struct_synth_th{3,5}.area,'ob'); hold on
plot(area_synth_v3_mem5,SSI_struct_synth_th_noRUSICE.area,'xk');
xlabel('Affected Area, reported by WISC [km^2]'); ylabel('Affected Area, calculated by CLIMADA [km^2]'); 

figure; plot(SSI_synth_v3_mem5,SSI_synth_th{3,5},'or'); hold on
plot(SSI_synth_v3_mem5,SSI_synth_th_noRUSICE,'xk'); 
xlabel('SSI, reported by WISC '); ylabel('SSI, calculated by CLIMADA '); 

figure; plot(gust_synth_v3_mem5, SSI_struct_synth_th{3,5}.gustspeed,'og'); hold on
plot(gust_synth_v3_mem5,SSI_struct_synth_th_noRUSICE.gustspeed,'xk');
xlabel('Mean gustspeed, reported by WISC [m/s]'); ylabel('Affected Area, calculated by CLIMADA [m/s]'); 

% => everything better but still not very satisfying
sum(hazard_wisc_synth_no_russia_iceland.on_land)
sum(hazard_wisc_synth.on_land)
% hypothesis 3: maybe both? russia, iceland and unit area
[SSI_synth_th_both, SSI_struct_synth_th_both] = climada_hazard_SSI_calc(hazard_wisc_synth_no_russia_iceland,25,ones(size(hazard_wisc_synth_no_russia_iceland.area_km2))*11.69);

figure; plot(SSI_synth_v3_mem5,SSI_synth_th{3,5},'or'); hold on
plot(SSI_synth_v3_mem5,SSI_synth_th_both,'xk'); 
% => the best fit but still not a perfect match (in single cases factor
% two) 
% next step look at event (geographical extent of intensities, climada_hazard_plot) with biggest difference for new ideas where the
% difference comes from...
% save('20181116_Wisc_new_synth_case.mat','-v7.3')

test_event_id = find((SSI_synth_v3_mem5 < 0.5*10^10) & ( SSI_synth_th_both' > 10^10)); %select the event with ssi difference of factor 2, difference between "wisc reported ssi" and calculated ssi
figure; climada_hazard_plot(hazard_wisc_synth,test_event_id);
% => this specific high intensities in northern africa and greenland
% hypothesis 4: with the exclusion of greenland and northern africa in the calculation of the ssi the fit
% with "wisc reported ssi should be minimal
hazard_wisc_synth_no_non_europe = hazard_wisc_synth; 
[country_name,~,~,incountry_russia] = climada_country_name('Russia','',hazard_wisc_synth_no_non_europe.lon,hazard_wisc_synth_no_non_europe.lat);
[country_name,~,~,incountry_iceland] = climada_country_name('Iceland','',hazard_wisc_synth_no_non_europe.lon,hazard_wisc_synth_no_non_europe.lat);
[country_name,~,~,incountry_greenland] = climada_country_name('Greenland','',hazard_wisc_synth_no_non_europe.lon,hazard_wisc_synth_no_non_europe.lat);
[country_name,~,~,incountry_morocco] = climada_country_name('Morocco','',hazard_wisc_synth_no_non_europe.lon,hazard_wisc_synth_no_non_europe.lat);
[country_name,~,~,incountry_tunisia] = climada_country_name('Tunisia','',hazard_wisc_synth_no_non_europe.lon,hazard_wisc_synth_no_non_europe.lat);
[country_name,~,~,incountry_algeria] = climada_country_name('Algeria','',hazard_wisc_synth_no_non_europe.lon,hazard_wisc_synth_no_non_europe.lat);
hazard_wisc_synth_no_non_europe.on_land(incountry_russia | incountry_iceland |...
    incountry_greenland | incountry_morocco | incountry_tunisia |...
    incountry_algeria) = false;
[SSI_synth_th_unitarea_noneurope, SSI_struct_synth_th_unitarea_noneurope] = climada_hazard_SSI_calc(hazard_wisc_synth_no_non_europe,25,ones(size(hazard_wisc_synth_no_non_europe.area_km2))*11.69);
% plot correlations
figure; plot(SSI_synth_v3_mem5,SSI_synth_th{3,5},'or'); hold on
plot(SSI_synth_v3_mem5,SSI_synth_th_unitarea_noneurope,'xk'); 
% => still not as perfect as hoped, specially because the order would not
% be identical with both ssi
% hypothesis 4: non europe excluded but regular area calculation
[SSI_synth_th_noneurope, SSI_struct_synth_th_noneurope] = climada_hazard_SSI_calc(hazard_wisc_synth_no_non_europe,25,hazard_wisc_synth_no_non_europe.area_km2);
% plot correlations
figure; plot(SSI_synth_v3_mem5,SSI_synth_th{3,5},'or'); hold on
plot(SSI_synth_v3_mem5,SSI_synth_th_noneurope,'xk'); 
% => fit very good, especially order of big events. still a few exceptions
% look at one event with large difference
test_event_id2 = find((SSI_synth_v3_mem5 < 2*10^10) & ( SSI_synth_th_noneurope' > 2*10^10)); %select the event with ssi difference of factor 2, difference between "wisc reported ssi" and calculated ssi
figure; climada_hazard_plot(hazard_wisc_synth,test_event_id2);
% => high intensity in syria, -> exclude syria additionally to other countries
[country_name,~,~,incountry_syria] = climada_country_name('Syria','',hazard_wisc_synth_no_non_europe.lon,hazard_wisc_synth_no_non_europe.lat);
hazard_wisc_synth_no_non_europe.on_land(incountry_syria) = false;

[SSI_synth_th_noneurope, SSI_struct_synth_th_noneurope] = climada_hazard_SSI_calc(hazard_wisc_synth_no_non_europe,25,hazard_wisc_synth_no_non_europe.area_km2);
% plot correlations
figure; plot(SSI_synth_v3_mem5,SSI_synth_th{3,5},'or'); hold on
plot(SSI_synth_v3_mem5,SSI_synth_th_noneurope,'xk'); 
xlabel('SSI, reported by WISC '); ylabel('SSI, calculated by CLIMADA '); 
title({'Exclusion of non-europe countries/landmasses';'like Russia, Greenland, Northern Africa and Syria'})
%=> reasonable fit and match achieved
% check area and mean gust
figure; plot(area_synth_v3_mem5, SSI_struct_synth_th{3,5}.area,'ob'); hold on
plot(area_synth_v3_mem5,SSI_struct_synth_th_noneurope.area,'xk');
xlabel('Affected Area, reported by WISC [km^2]'); ylabel('Affected Area, calculated by CLIMADA [km^2]'); 
title({'Exclusion of non-europe countries/landmasses';'like Russia, Greenland, Northern Africa and Syria'})
% area matches well
figure; plot(gust_synth_v3_mem5, SSI_struct_synth_th{3,5}.gustspeed,'og'); hold on
plot(gust_synth_v3_mem5,SSI_struct_synth_th_noneurope.gustspeed,'xk');
xlabel('Mean gustspeed, reported by WISC [m/s]'); ylabel('Affected Area, calculated by CLIMADA [m/s]'); 
title({'Exclusion of non-europe countries/landmasses';'like Russia, Greenland, Northern Africa and Syria'})
% mean gust does not match as well, why?
figure; plot(gust_synth_v3_mem5, SSI_struct_synth_th{3,5}.gustspeed,'og'); hold on
plot(gust_synth_v3_mem5,SSI_struct_synth_th_unitarea_noneurope.gustspeed,'xk');
% check davids method as cross comparison:
[SSI_synth_db_noneurope] = climada_hazard_ssi(hazard_wisc_synth_no_non_europe,0,25);
figure; plot(SSI_synth_v3_mem5,SSI_synth_db_noneurope,'xk'); 
% => does fit as well
test_event_id3 = find((SSI_synth_v3_mem5 < 1.65*10^10) & ( SSI_synth_db_noneurope' > 1.78*10^4)); %select the event with ssi difference of factor 2, difference between "wisc reported ssi" and calculated ssi
figure; climada_hazard_plot(hazard_wisc_synth,test_event_id3);
figure; plot(SSI_synth_th_noneurope,SSI_synth_db_noneurope,'xk'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% => Conclusion:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('The difference in SSI calculations between WISC and CLIMADA originates from the inclusion/exclusion of non-europe countries/landmasses like Russia, Greenland, Northern Africa and Syria.\n') 


