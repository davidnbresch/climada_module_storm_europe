%function climatewise_run
% climada climatewise
% MODULE:
%   storm europe
% NAME:
%   climatewise_run
% PURPOSE:
%   batch job to analise Climate Wise data - see climatewise_core for
%   the core batch job which is called by climatewise_run
%
%   Process the GBR_*.xlsx exposure data in one batch and the commercial 
%   ones ??.xlsx in a second one (see exposure_files in PARAMETERS below)
%
%   for GBR_barclays_sector_agg.xlsx etc:
%   process_number_of_properties=1;
%   Percentage_Of_Value_Flag=1;
%
%   for ??.xlsx  etc:
%   process_number_of_properties=0;
%   Percentage_Of_Value_Flag=0;
%
%   Input are the bank's assets, stored in single Excel files with one
%   header row, followed by data (all numeric except for postcode_sector,
%   where AB10 1 is the entry, hence this space not to be confused with
%   space as a separator, as it might appear in the raw text here, but all
%   fine within Excel):
%       postcode_sector	latitude	longitude       number of properties    replacement_value_gbp
%       AB10 1          57.14687255	-2.106558026	17                      3447974
%       AB10 6          57.13707772	-2.122704986	106                     21499132
%       ...
%
%   Please note that for speedup, the hazard once loaded is kept, hence run
%   clear hazard in case you switch to another hazard set. Read the
%   PARAMETERS section carefully anyway (an expert code, not for beginners)
%   Please familiarize yurself with CLIMADA before running this code.
%
% CALLING SEQUENCE:
%   climatewise_run
% EXAMPLE:
%   climatewise_run
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to stdout and figures and to a folder ClimateWise within the CLIMADA
%   results folder
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20170803, initial, separated from climatewise_core
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
%
% most parameters are set below, here only the general ones
%
% the reference return period we express the results for
reference_return_period=100;
%
Intensity_threshold_ms_WS=35; % intensity threshold for affected in m/s for WS
Intensity_threshold_ms_TC=55; % intensity threshold for affected in m/s for TC
%
exposure_folder=[climada_global.data_dir filesep 'ClimateWise'];


% first GBR European winter storm only
% ====================================
    
% first for risk today
% --------------------
clear hazard % to be on the safe side
fprintf('\n\n*** GBR_*.xlsx risk today ***\n\n');

climate_frequency_screw=0; % indicates no climate change
exposure_files={
    'GBR_barclays_sector_agg.xlsx' % GBR_ to indicate GBR only (speedup, must be listed first)
    'GBR_hsbc_sector_agg.xlsx'
    'GBR_lloyds_sector_agg.xlsx'
    'GBR_nationwide_sector_agg.xlsx'
    'GBR_rbs_sector_agg.xlsx'
    'GBR_santander_sector_agg.xlsx'
    'GBR_yorkshire_sector_agg.xlsx'
    };
process_number_of_properties=1;
Percentage_Of_Value_Flag=1;
climatewise_core % call the core function

% store results
DFC=climada_EDS2DFC(EDS); % calculate all damage excedance frequency curves
reference_pos=find(DFC(1).return_period==reference_return_period);

for DFC_i=1:length(DFC)
    GBR_result_struct(DFC_i).annotation_name=DFC(DFC_i).annotation_name;
    GBR_result_struct(DFC_i).risk_today=DFC(DFC_i).damage(reference_pos)/DFC(DFC_i).Value*100;
end

% second for risk climate change
% ------------------------------
clear hazard % to be on the safe side
fprintf('\n\n*** GBR_*.xlsx risk climate change ***\n\n');

climate_frequency_screw=1; % indicates climate change
climate_intensity_a=1.0;climate_intensity_b=0.2; % intensity=a*intensity+b

climatewise_core % call the core function

% store results
DFC=climada_EDS2DFC(EDS); % calculate all damage excedance frequency curves
reference_pos=find(DFC(1).return_period==reference_return_period);

for DFC_i=1:length(DFC)
    if strcmpi(DFC(DFC_i).annotation_name,GBR_result_struct(DFC_i).annotation_name)
        GBR_result_struct(DFC_i).risk_climate=DFC(DFC_i).damage(reference_pos)/DFC(DFC_i).Value*100;
        GBR_result_struct(DFC_i).risk_delta=(GBR_result_struct(DFC_i).risk_climate/GBR_result_struct(DFC_i).risk_today-1)*100;
    else
        fprintf('ERROR: exposure_files mismatch, aborted\n');
        return
    end
end

fprintf('\nexposure: \t\t today \t\t climate \t\t delta\n ');
for DFC_i=1:length(DFC)
    fprintf('%s: \t\t %2.0f%% \t\t %2.0f%% \t\t %2.0f%%\n ',...
        GBR_result_struct(DFC_i).annotation_name,...
        GBR_result_struct(DFC_i).risk_today,...
        GBR_result_struct(DFC_i).risk_climate,...
        GBR_result_struct(DFC_i).risk_delta);
end



% second commercial exposures (WS and TC)
% =======================================
    
% first for risk today
% --------------------
clear hazard % to be on the safe side
fprintf('\n\n*** ??.xlsx risk today ***\n\n');

climate_frequency_screw=0; % indicates no climate change
exposure_files={
    '02.xlsx' % follow the global exposures
    '03.xlsx'
    '05.xlsx'
    '08.xlsx'
    '12.xlsx'
    '13.xlsx'
    '14.xlsx'
    '15.xlsx'
    '16.xlsx'
    '23.xlsx'
    '24.xlsx'
    '25.xlsx'
    };
process_number_of_properties=0;
Percentage_Of_Value_Flag=0;
climatewise_core % call the core function

% store results
DFC_WS=climada_EDS2DFC(EDS_WS); % calculate all damage excedance frequency curves
reference_pos=find(DFC_WS(1).return_period==reference_return_period);

for DFC_i=1:length(DFC_WS)
    XX_WS_result_struct(DFC_i).annotation_name=DFC_WS(DFC_i).annotation_name;
    XX_WS_result_struct(DFC_i).risk_today=DFC_WS(DFC_i).damage(reference_pos)/DFC_WS(DFC_i).Value*10000;
    XX_WS_result_struct(DFC_i).risk_climate=0;
    XX_WS_result_struct(DFC_i).risk_delta=NaN;
end

DFC_TC=climada_EDS2DFC(EDS_TC); % calculate all damage excedance frequency curves
reference_pos=find(DFC_TC(1).return_period==reference_return_period);

for DFC_i=1:length(DFC_TC)
    XX_TC_result_struct(DFC_i).annotation_name=DFC_TC(DFC_i).annotation_name;
    XX_TC_result_struct(DFC_i).risk_today=DFC_TC(DFC_i).damage(reference_pos)/DFC_TC(DFC_i).Value*10000;
    XX_TC_result_struct(DFC_i).risk_climate=0;
    XX_TC_result_struct(DFC_i).risk_delta=NaN;
end

fprintf('\nWS exposure \t\t today \t\t climate \t\t delta\n ');
for DFC_i=1:length(XX_WS_result_struct)
    fprintf('%s: \t\t %2.0f%%oo \t\t %2.0f%%oo \t\t %2.0f%%\n ',...
        XX_WS_result_struct(DFC_i).annotation_name,...
        XX_WS_result_struct(DFC_i).risk_today,...
        XX_WS_result_struct(DFC_i).risk_climate,...
        XX_WS_result_struct(DFC_i).risk_delta);
end

fprintf('\n TC exposure \t\t today \t\t climate \t\t delta\n ');
for DFC_i=1:length(XX_TC_result_struct)
    fprintf('%s: \t\t %2.0f%%oo \t\t %2.0f%%oo \t\t %2.0f%%\n ',...
        XX_TC_result_struct(DFC_i).annotation_name,...
        XX_TC_result_struct(DFC_i).risk_today,...
        XX_TC_result_struct(DFC_i).risk_climate,...
        XX_TC_result_struct(DFC_i).risk_delta);
end

