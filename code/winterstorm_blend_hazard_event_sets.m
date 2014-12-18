function hazard=winterstorm_blend_hazard_event_sets(hazard_save_file,frequency_screw,intensity_screw)
% climada
% MODULE:
%   ws_europe
% NAME:
%   winterstorm_blend_hazard_event_sets
% PURPOSE:
%   Blend the four WS hazard events sets into one, see PARAMETERS in code
%
%   Hazard sets from:
%   Schwierz, C., P. Koellner-Heck, E. Zenklusen Mutter, D. N. Bresch,
%   P.-L.Vidale, M. Wild, C., and Sch?r, 2010: Modelling European winter
%   wind storm losses in current and future climate. Climatic Change (2010)
%   101:485?514, doi: 10.1007/s10584-009-9712-1.
%
%   See also winterstorm_compare and winterstorm_compare_severity
% CALLING SEQUENCE:
%   hazard=winterstorm_blend_hazard_event_sets(hazard_save_file,frequency_screw)
% EXAMPLE:
%   hazard=winterstorm_blend_hazard_event_sets('WS_Europe_blend.mat',0.9)
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   hazard_save_file: the filename (w/o path) of the blended output file
%       It is always stored in the same folder as the files it is blended
%       from, e.g. only a filename, such as 'WS_Europe_blend.mat'.
%   frequency_screw: an EXPERIMENTAL multiplier for the frequency, just
%       multiplies all single event frequencies. Default=1 (obviously)
%       (Can be justified to adjust, due to the fact that we blend hazard
%       sets, an adjustment accounts for especially smaller scale events) 
%   intensity_screw: scale intensity with, default=1
% OUTPUTS:
%   hazard: the blended hazard event set, usually named WS_Europe.mat, see
%       optional input parameter hazard_save_file
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141201, initial
% David N. Bresch, david.bresch@gmail.com, 20141206, frequency_screw added
% David N. Bresch, david.bresch@gmail.com, 20141213, intensity_screw added
%-

hazard=[]; % init output
hazard_blend=[];

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('hazard_save_file','var'),hazard_save_file='';end
if ~exist('frequency_screw','var'),frequency_screw=1;end % default=1
if ~exist('intensity_screw','var'),intensity_screw=1;end % default=1

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
hazard_set_folder=[module_data_dir filesep 'hazards'];
hazard_set_files={'WS_ECHAM_CTL','WS_ETHC_CTL','WS_GKSS_CTL','WS_ERA40'};
if isempty(hazard_save_file)
    hazard_save_file=[hazard_set_folder filesep 'WS_Europe.mat'];
else
    hazard_save_file=[hazard_set_folder filesep hazard_save_file];
end

for file_i=1:length(hazard_set_files)
    hazard_set_file=[hazard_set_folder filesep hazard_set_files{file_i} '.mat'];
    if exist(hazard_set_file,'file')
        
        fprintf('adding %s\n',hazard_set_files{file_i});
        
        load(hazard_set_file);
        
        if isempty(hazard_blend)
            hazard_blend=hazard;
            hazard.comment=['WSEU blended ' hazard_set_files{file_i}];
            hazard.date=datestr(now);
            hazard_blend.hazard_count=1;
        else
            hazard_blend.intensity  =[hazard_blend.intensity; hazard.intensity];
            hazard_blend.comment    =[hazard_blend.comment ' ' hazard_set_files{file_i}];
            hazard_blend.event_count=hazard_blend.event_count+hazard.event_count;
            hazard_blend.orig_event_flag=[hazard_blend.orig_event_flag hazard.orig_event_flag];
            hazard_blend.frequency  =[hazard_blend.frequency hazard.frequency];
            hazard_blend.event_ID   =1:hazard_blend.event_count;
            hazard_blend.orig_years =max(hazard_blend.orig_years,hazard.orig_years);
            hazard_blend.filename   =sprintf('blended, see %s',mfilename);
            hazard_blend.hazard_count=hazard_blend.hazard_count+1;
        end % file_i==1
    end % exist(hazard_set_file)
end % file

hazard_blend.frequency=hazard_blend.frequency/hazard_blend.hazard_count*frequency_screw;
hazard_blend.frequency_screw_applied=frequency_screw;

hazard=hazard_blend;
hazard.intensity(hazard.intensity<15)=0;
% apply intensity_screw after cutoff, in order to allow for 'turn back' without loosing events
hazard.intensity=sparse(hazard.intensity)*intensity_screw; % sparsify again, to be safe
hazard.intensity_screw_applied=intensity_screw;
hazard.matrix_density = nnz(hazard.intensity)/numel(hazard.intensity);

save(hazard_save_file,'hazard');
climada_hazard_cleanup(hazard_save_file);

return
