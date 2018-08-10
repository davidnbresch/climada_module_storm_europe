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
% David N. Bresch, david.bresch@gmail.com, 20170810, all automatic for WS and TC, today and climate change
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
wsgsmax_diff_file = [climada_global.data_dir filesep 'ClimateWise' filesep 'wsgsmax_diff.nc'];
% since the times on the netCDF file wsgsmax_diff.nc are from a
% non-standard starting point, we hard-wire them:
WS_climate_scenario_times=[2015 20125 2035 2045]; % reference 1960-1990
WS_time_i=4; % for timestep 1 on wsgsmax_diff.nc
%
Intensity_threshold_ms_WS=35; % intensity threshold for affected in m/s for WS
Intensity_threshold_ms_TC=55; % intensity threshold for affected in m/s for TC
%
exposure_folder=[climada_global.data_dir filesep 'ClimateWise'];


% first GBR European winter storm only
% ====================================

exposure_files={
    'GBR_barclays_sector_agg.xlsx' % GBR_ to indicate GBR only (speedup, must be listed first)
    'GBR_hsbc_sector_agg.xlsx'
    'GBR_lloyds_sector_agg.xlsx'
    'GBR_nationwide_sector_agg.xlsx'
    'GBR_rbs_sector_agg.xlsx'
    'GBR_santander_sector_agg.xlsx'
    'GBR_yorkshire_sector_agg.xlsx' % for tests, use only first and last
    };

process_number_of_properties=1;
Percentage_Of_Value_Flag=1;

% first for risk today
% --------------------
clear hazard % to be on the safe side
WS_hazard_CC_ext=''; % no climate change
fprintf('\n\n*** GBR_*.xlsx risk today ***\n\n');

climatewise_core % call the core function

% store results for risk today
DFC=climada_EDS2DFC(EDS); % calculate all damage excedance frequency curves
reference_pos = DFC(1).return_period==reference_return_period;

for DFC_i=1:length(DFC)
    GBR_result_struct(DFC_i).annotation_name=DFC(DFC_i).annotation_name;
    if Percentage_Of_Value_Flag
        GBR_result_struct(DFC_i).risk_today=DFC(DFC_i).damage(reference_pos)/DFC(DFC_i).Value*100;
    else
        GBR_result_struct(DFC_i).risk_today=DFC(DFC_i).damage(reference_pos);
    end
end

% second for risk climate change
% ------------------------------
clear hazard % to be on the safe side
WS_hazard_CC_ext = sprintf('_CC%4.4i',WS_climate_scenario_times(WS_time_i)); % make sure WS_time_i matches
fprintf('\n\n*** GBR_*.xlsx risk climate change ***\n\n');

climatewise_core % call the core function

% store results for climate change
DFC=climada_EDS2DFC(EDS); % calculate all damage excedance frequency curves
reference_pos = DFC(1).return_period==reference_return_period;

for DFC_i=1:length(DFC)
    if strcmpi(DFC(DFC_i).annotation_name,GBR_result_struct(DFC_i).annotation_name)
        if Percentage_Of_Value_Flag
            GBR_result_struct(DFC_i).risk_climate=DFC(DFC_i).damage(reference_pos)/DFC(DFC_i).Value*100;
        else
            GBR_result_struct(DFC_i).risk_climate=DFC(DFC_i).damage(reference_pos);
        end
        GBR_result_struct(DFC_i).risk_delta=(GBR_result_struct(DFC_i).risk_climate/GBR_result_struct(DFC_i).risk_today-1)*100;
    else
        fprintf('ERROR: exposure_files mismatch, aborted\n');
        return
    end
end

% figure what we display
if     log10(EDS(1).currency_unit)==0
    unit_str=EDS(1).Value_unit;
elseif log10(EDS(1).currency_unit)==6
    unit_str=[EDS(1).Value_unit ' mio'];
elseif log10(EDS(1).currency_unit)==9
    unit_str=[EDS(1).Value_unit ' bn'];
end
if process_number_of_properties==1
    title_str='number of properties affected';
else
    title_str='values';
end
if Percentage_Of_Value_Flag==1
    title_str=[title_str ' in percent'];
    fmt_str='%s: \t\t %2.0f%% \t\t %2.0f%% \t\t %2.1f%%\n';
else
    if ~process_number_of_properties
        title_str=[title_str ' in ' unit_str];
    end
    fmt_str='%s: \t\t %2.3f \t\t %2.3f \t\t %2.1f%%\n';
end
fprintf('%s\n',title_str);

fprintf('\nexposure: \t\t today \t\t climate \t\t delta\n');
for DFC_i=1:length(DFC)
    fprintf(fmt_str,...
        GBR_result_struct(DFC_i).annotation_name,...
        GBR_result_struct(DFC_i).risk_today,...
        GBR_result_struct(DFC_i).risk_climate,...
        GBR_result_struct(DFC_i).risk_delta);
end


% number of properties affected in percent
% exposure: 		 today 		 climate     delta
% barclays: 		 43% 		 45% 		 4.9%
% hsbc:              40% 		 42% 		 5.8%
% lloyds:            39% 		 41% 		 3.9%
% nationwide:        40% 		 42% 		 5.5%
% rbs:               39% 		 41% 		 6.1%
% santander: 		 40% 		 42% 		 5.1%
% yorkshire: 		 38% 		 39% 		 4.3%
% combined: 		 38% 		 41% 		 5.6%


% second commercial exposures (WS and TC)
% =======================================

exposure_files={
    '02.xlsx' % follow the global exposures
    '03.xlsx'
    '05.xlsx'
    '08.xlsx' % use this for tests
    '12.xlsx'
    '13.xlsx'
    '14.xlsx'
    '15.xlsx'
    '16.xlsx'
    '23.xlsx'
    '24.xlsx'
    '25.xlsx' % and this for tests, too
    };

process_number_of_properties=0;
Percentage_Of_Value_Flag=0;

% first for risk today
% --------------------
clear hazard % to be on the safe side
WS_hazard_CC_ext=''; % no climate change
TC_hazard_CC_ext=''; % no climate change
fprintf('\n\n*** ??.xlsx risk today ***\n\n');

climatewise_core % call the core function

% store results for risk today
DFC_WS=climada_EDS2DFC(EDS_WS); % calculate all damage excedance frequency curves
reference_pos= DFC_WS(1).return_period==reference_return_period;

for DFC_i=1:length(DFC_WS)
    XX_WS_result_struct(DFC_i).annotation_name=DFC_WS(DFC_i).annotation_name;
    XX_WS_result_struct(DFC_i).risk_today=DFC_WS(DFC_i).damage(reference_pos);
end

DFC_TC=climada_EDS2DFC(EDS_TC); % calculate all damage excedance frequency curves
reference_pos= DFC_TC(1).return_period==reference_return_period;

for DFC_i=1:length(DFC_TC)
    XX_TC_result_struct(DFC_i).annotation_name=DFC_TC(DFC_i).annotation_name;
    XX_TC_result_struct(DFC_i).risk_today=DFC_TC(DFC_i).damage(reference_pos);
end

% second for climate change
% -------------------------
clear hazard % to be on the safe side
WS_hazard_CC_ext='_CC2015';     % climate change
TC_hazard_CC_ext='_rcp45_2030'; % climate change
fprintf('\n\n*** ??.xlsx climate change ***\n\n');

climatewise_core % call the core function

% store results for climate change
DFC_WS=climada_EDS2DFC(EDS_WS); % calculate all damage excedance frequency curves
reference_pos = DFC_WS(1).return_period==reference_return_period;

for DFC_i=1:length(DFC_WS)
    if strcmpi(DFC_WS(DFC_i).annotation_name,XX_WS_result_struct(DFC_i).annotation_name)
        XX_WS_result_struct(DFC_i).risk_climate=DFC_WS(DFC_i).damage(reference_pos);
        XX_WS_result_struct(DFC_i).risk_delta=(XX_WS_result_struct(DFC_i).risk_climate/XX_WS_result_struct(DFC_i).risk_today-1)*100;
    else
        fprintf('ERROR: exposure_files mismatch, aborted\n');
        return
    end
end % DFC_i for WS

DFC_TC=climada_EDS2DFC(EDS_TC); % calculate all damage excedance frequency curves
reference_pos = DFC_TC(1).return_period==reference_return_period;

for DFC_i=1:length(DFC_TC)
    if strcmpi(DFC_TC(DFC_i).annotation_name,XX_TC_result_struct(DFC_i).annotation_name)
        XX_TC_result_struct(DFC_i).risk_climate=DFC_TC(DFC_i).damage(reference_pos);
        XX_TC_result_struct(DFC_i).risk_delta=(XX_TC_result_struct(DFC_i).risk_climate/XX_TC_result_struct(DFC_i).risk_today-1)*100;
    else
        fprintf('ERROR: exposure_files mismatch, aborted\n');
        return
    end
end % DFC_i for TC

% figure what we display
if     log10(EDS_WS(1).currency_unit)==0
    unit_str=EDS_WS(1).Value_unit;
elseif log10(EDS_WS(1).currency_unit)==6
    unit_str=[EDS_WS(1).Value_unit ' mio'];
elseif log10(EDS_WS(1).currency_unit)==9
    unit_str=[EDS_WS(1).Value_unit ' bn'];
end

fprintf('\nrisk in %s\n',unit_str);
fprintf('WS exposure \t today \t\t climate \t\t delta\n');
for DFC_i=1:length(XX_WS_result_struct)
    fprintf('%s: \t\t %2.3f \t\t %2.3f \t\t %2.1f%%\n ',...
        XX_WS_result_struct(DFC_i).annotation_name,...
        XX_WS_result_struct(DFC_i).risk_today,...
        XX_WS_result_struct(DFC_i).risk_climate,...
        XX_WS_result_struct(DFC_i).risk_delta);
end

fprintf('TC exposure \t today \t\t climate \t\t delta\n');
for DFC_i=1:length(XX_TC_result_struct)
    fprintf('%s: \t\t %2.3f \t\t %2.3f \t\t %2.1f%%\n ',...
        XX_TC_result_struct(DFC_i).annotation_name,...
        XX_TC_result_struct(DFC_i).risk_today,...
        XX_TC_result_struct(DFC_i).risk_climate,...
        XX_TC_result_struct(DFC_i).risk_delta);
end

% risk in GBP mio
% WS exposure today 	 climate 	 delta
%  02: 		 0.241 		 0.227 		 -5.6%
%  03: 		 1.907 		 2.035 		 6.7%
%  05: 		 0.191 		 0.186 		 -2.9%
%  08: 		 0.328 		 0.344 		 4.7%
%  12: 		 0.512 		 0.524 		 2.2%
%  13: 		 0.017 		 0.018 		 7.3%
%  15: 		 0.008 		 0.008 		 6.7%
%  16: 		 0.056 		 0.056 		 -0.3%
%  23: 		 0.090 		 0.096 		 6.4%
%  24: 		 0.000 		 0.000 		 5.9%
%  25: 		 0.536 		 0.588 		 9.8%
%  combined: 3.602 		 3.801 		 5.5%
%
% TC exposure today      climate     delta
%  08: 		 0.355 		 0.431 		 21.5%
%  12: 		 0.914 		 1.155 		 26.4%
%  14: 		 0.451 		 0.559 		 24.0%
%  15: 		 0.002 		 0.002 		 26.6%
%  25: 		 0.475 		 0.578 		 21.8%
%  combined: 1.566 		 1.890 		 20.7%
