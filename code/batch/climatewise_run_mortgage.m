%function climatewise_run_mortgage
% climada climatewise
% MODULE:
%   storm europe
% NAME:
%   climatewise_run_mortgage
% PURPOSE:
%   copy of climatewise_run, just for the first part, see there. Best run
%   in Octave. See climatewise_run first.
%
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
% David N. Bresch, david.bresch@gmail.com, 20180803, initial, separated from climatewise_core
% Mark Westcott, mark.westcott@vivideconomics.com, 20180901, spearated from climatewise_run
%-

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

%% Climate change scenarios

climate_scenario_times=[2015 2025 2035 2045 2055]; % reference 1960-1990


%# Set up summary output file

summary_filename = [output_dir filesep 'mortgage-summary.csv'];
fid = fopen(summary_filename,'w'); 
fprintf(fid,'portfolio,scenario,peril,year,ED_value,risk_today,value_column,currency_units,adaptation\n');
fclose(fid);

%%%% Cells for summary results, later written to a CSV
portfolio_column = {};
scenario_column = {};
peril_column = {};
year_column = {};
ED_value_column = {};
value_units_column = {};
currency_units_column = {};
risk_today_column = {};
value_column = {};
adaptation_column = {};

% Clear out (to be safe)
clear EDS_WS;
clear EDS_local;



% first for baseline today
% --------------------

fprintf('\n\n*** GBR_*.xlsx baseline risk ***\n\n');


clear hazard % to be on the safe side

WS_hazard_CC_ext=''; % no climate change
scenario_name = 'baseline';
year = 1975;

adaptation = ''; % no adaptation
climatewise_core % call the core function
climatewise_write_eds_to_csv %write out ED curves to CSV and populate cells with summary results



% RCP 4.5 from 2015 to 2055
% ------------------------------


wsgsmax_diff_file = [climada_global.data_dir filesep 'ClimateWise' filesep 'rcp45-wsgsmax-99pctl-diff.nc'];
scenario_name = 'rcp45';

for time_i=1:length(climate_scenario_times)
  clear hazard % to be on the safe side
  year = climate_scenario_times(time_i);
  WS_time_i = time_i; % Which time period to read from wsgsmax_diff_file
  WS_hazard_CC_ext = sprintf('_rcp45_%4.4i',year);
  
  fprintf('\n\n*** GBR_*.xlsx RCP 4.5 ***\n\n');

  %# Without adaptation 
  adaptation = '';
  climatewise_core % call the core function
  climatewise_write_eds_to_csv %write out ED curves to CSV and populate cells with summary results
  
%   # With adaptation
%   #adaptation = '-adaptation';
%   #climatewise_core % call the core function
%   #climatewise_write_eds_to_csv

end

return;




% RCP 8.5 from 2015 to 2055
% ------------------------------

wsgsmax_diff_file = [climada_global.data_dir filesep 'ClimateWise' filesep 'rcp85-wsgsmax-99pctl-diff.nc'];

scenario_name = 'rcp85';

for time_i=1:length(climate_scenario_times)
  clear hazard % to be on the safe side
  year = climate_scenario_times(time_i);
  WS_time_i = time_i; % Which time period to read from wsgsmax_diff_file
  WS_hazard_CC_ext = sprintf('_rcp85_%4.4i',year);

  fprintf('\n\n*** GBR_*.xlsx risk climate change RCP 8.5 ***\n\n');

  climatewise_core % call the core function
  
  %# Without adaptation 
  adaptation = '';
  climatewise_core % call the core function
  climatewise_write_eds_to_csv
  
%   # With adaptation
%   #adaptation = '-adaptation';
%   #climatewise_core % call the core function
%   #climatewise_write_eds_to_csv

end


%#Write summary

fid = fopen(summary_filename,'a'); 
for i=1:length(portfolio_column)

  formatted_string = sprintf ("%s,%s,%s,%d,%d,%d,%s,%d,%s\n",
        portfolio_column{i}, scenario_column{i}, peril_column{i}, year_column{i}, ED_value_column{i},  risk_today_column{i} ,value_column{i},  currency_units_column{i}, adaptation_column{i});
  fprintf(fid, formatted_string);
end
fclose(fid);