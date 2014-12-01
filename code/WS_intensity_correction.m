function WS_intensity_correction(save_files)
% climada
% NAME:
%   WS_intensity_correction
% PURPOSE:
%   correct WS intenisty as documented in the paper:
%
%   Schwierz, C., P. K?llner-Heck, E. Zenklusen Mutter, D. N. Bresch,
%   P.-L.Vidale, M. Wild, C., and Sch?r, 2010: Modelling European winter
%   wind storm losses in current and future climate. Climatic Change (2010)
%   101:485?514, doi: 10.1007/s10584-009-9712-1.
%
%   If applied a second time, the correction is reversed
%
%   WARNING: a truly expert level code, to be used with utmost caution
%
%   see also winterstorm_compare and winterstorm_compare_severity
% CALLING SEQUENCE:
%   WS_intensity_correction(param1);
% EXAMPLE:
%   WS_intensity_correction
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   save_files: whether files are saves(=1) or only inspected (=0, default)
% OUTPUTS:
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141129
%-

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout
if ~exist('save_files','var'),save_files=0;end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
hazard_set_folder=[module_data_dir filesep 'hazards'];
hazard_set_files= {'WS_ECHAM_CTL','WS_ECHAM_A2','WS_ETHC_CTL','WS_ETHC_A2','WS_GKSS_CTL','WS_GKSS_A2'};
hazard_set_factors=[1.22           1.22          1.17          1.17         1.09          1.09];
intensity_comment{1}='intensity multiplied by factor 1.22 according to Schwierz et al, 2010, doi: 10.1007/s10584-009-9712-1';
intensity_comment{2}='intensity multiplied by factor 1.22 according to Schwierz et al, 2010, doi: 10.1007/s10584-009-9712-1';
intensity_comment{3}='intensity multiplied by factor 1.17 (HC-CHRM) according to Schwierz et al, 2010, doi: 10.1007/s10584-009-9712-1';
intensity_comment{4}='intensity multiplied by factor 1.17 (HC-CHRM) according to Schwierz et al, 2010, doi: 10.1007/s10584-009-9712-1';
intensity_comment{5}='intensity multiplied by factor 1.09 (HC-CLM) according to Schwierz et al, 2010, doi: 10.1007/s10584-009-9712-1';
intensity_comment{6}='intensity multiplied by factor 1.09 (HC-CLM) according to Schwierz et al, 2010, doi: 10.1007/s10584-009-9712-1';

fprintf('check %i hazard event sets:\n',length(hazard_set_files)+1);

for hazard_i=1:length(hazard_set_files)
    hazard_set_file=[hazard_set_folder filesep hazard_set_files{hazard_i}];
    hazard_set_short=strrep(hazard_set_files{hazard_i},'WS_','');
    fprintf('%s:\n',hazard_set_short);
    load(hazard_set_file)
    
    if isfield(hazard,'intensity_comment')
        fprintf('on file: %s (%2.2f)\n',hazard.intensity_comment,hazard_set_factors(hazard_i));
        hazard.intensity=hazard.intensity/hazard_set_factors(hazard_i);
        hazard.intensity=sparse(hazard.intensity);
        hazard=rmfield(hazard,'intensity_comment');
        if save_files
            fprintf('intensity reset and saved (overwritten)\n');
            save(hazard_set_file,'hazard');
        end
    else
        hazard.intensity_comment=intensity_comment{hazard_i};
        fprintf('adding: %s (%2.2f)\n',hazard.intensity_comment,hazard_set_factors(hazard_i));
        hazard.intensity=hazard.intensity*hazard_set_factors(hazard_i);
        hazard.intensity=sparse(hazard.intensity);
        if save_files
            fprintf('intensity adjusted and saved (overwritten)\n');
            save(hazard_set_file,'hazard');
        end
    end

end % hazard_i


return
