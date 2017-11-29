function res=schwierz_etal_2010(return_period)
% climada schwierz 2010
% MODULE:
%   storm_europe
% NAME:
%   schwierz_etal_2010
% PURPOSE:
%   reproduce results of Schwierz et al, 2010. Plot a figure similar to
%   Fig. 10 in the paper.
%
%   Schwierz, C., K?llner-Heck, P., Zenklusen, E., Bresch, D. N., Vidale,
%   P.L., Wild, M., & Sch?r, C., 2010: Modelling European winter windstorm
%   losses in current and future climate, Climatic change, Vol. 101,
%   485-514. DOI 10.1007/s10584-009-9712-1.
%   http://www.iac.ethz.ch/doc/publications/Schwierz.pdf    
%
% CALLING SEQUENCE:
%   res=schwierz_etal_2010(return_period)
% EXAMPLE:
%   res=schwierz_etal_2010
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   return_period: the return period (e.g.=100) for which we show results,
%       default='AEL'
% OUTPUTS:
%   res: the results, empty if not successful
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20170505, initial
% David N. Bresch, david.bresch@gmail.com, 20170507, figure as in paper
% David N. Bresch, david.bresch@gmail.com, 20171129, entities generated
%-

res=[]; % init output

if ~exist('return_period','var'),return_period='AED';end

% global climada_global
% if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
climate_model={'ETHC','ECHAM','GKSS'};
country_ISO3={'DNK','DEU','SWE','BEL','NLD','FRA','GBR','CHE','NOR','IRL'};

n_countries=length(country_ISO3);
n_climate_models=length(climate_model);

% init
res.CTL=zeros(n_countries+1,n_climate_models);
res.A2 =res.CTL;

% construct the asset base
for country_i=1:n_countries
    entity=climada_entity_load(country_ISO3{country_i},-1); % try to load, silent mode
    if isempty(entity),entity=climada_entity_country(country_ISO3{country_i});end % construct it
end

% construct the combined entity
entity_EUR=climada_entity_load(country_ISO3{1});
for country_i=2:n_countries
    entity=climada_entity_load(country_ISO3{country_i});
    entity_EUR=climada_entity_combine(entity_EUR,entity);
end
% p.blue_ocean=1;climada_entity_plot(entity_EUR,1,p) % checkplot

for country_i=1:n_countries+1
    for climate_model_i=1:n_climate_models
        
        hazard_file_CTL=[module_data_dir filesep 'hazards' filesep 'WS_' climate_model{climate_model_i} '_CTL'];
        hazard_file_A2 =strrep(hazard_file_CTL,'_CTL','_A2');
        
        
        if country_i<=n_countries
            fprintf('processing %s for %s\n',country_ISO3{country_i},climate_model{climate_model_i});
            entity=climada_entity_load(country_ISO3{country_i});
        else
            fprintf('processing combined (EUR) for %s\n',climate_model{climate_model_i});
            entity=entity_EUR;
        end
        
        % to check with 'operational' WS Europe hazard set in climada
        % hazard_file_name=[entity.assets.admin0_ISO3 '_' strrep(entity.assets.admin0_name,' ','') '_eur_WS'];
        % hazard=climada_hazard_load(hazard_file_name);
        % EDS=climada_EDS_calc(entity,hazard);
        
        hazard_CTL=climada_hazard_load(hazard_file_CTL); % load CTL
        hazard_A2 =climada_hazard_load(hazard_file_A2);  % load A2
        
        entity_enc=climada_assets_encode(entity,hazard_CTL); % encode to low res
        EDS_CTL=climada_EDS_calc(entity_enc,hazard_CTL);
        EDS_A2 =climada_EDS_calc(entity_enc,hazard_A2);
        DFC_CTL=climada_EDS2DFC(EDS_CTL,return_period);
        DFC_A2 =climada_EDS2DFC(EDS_A2,return_period);
        
        res.CTL(country_i,climate_model_i)=DFC_CTL.damage;
        res.A2(country_i ,climate_model_i)=DFC_A2.damage;
        
    end % climate_model_i
    
end % country_i

% and now combined (EUR) first
res.CTL=circshift(res.CTL,1);
res.A2 =circshift(res.A2 ,1);

% calculate percentage difference
res.climate_signal=(res.A2-res.CTL)./res.CTL*100;

res.XTickLabel{1}='EUR';
res.XTickLabel(2:n_countries+1)=country_ISO3;

figure('Name','each model');
bar(res.climate_signal);set(gcf,'Color',[1 1 1])
set(gca,'XTickLabel',res.XTickLabel);
ylabel('Difference [% of CTL]');ylim([-50 150]);
title(climate_model)

% and the figure as close to the paper as possible (without fiddling ;-)
res.climate_signal_mean=mean(res.climate_signal,   2);
res.climate_signal_min = min(res.climate_signal,[],2)-res.climate_signal_mean;
res.climate_signal_max = max(res.climate_signal,[],2)-res.climate_signal_mean;

figure('Name','mean and extremes');
bar(res.climate_signal_mean,0.5,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0]);set(gcf,'Color',[1 1 1])
set(gca,'XTickLabel',res.XTickLabel);
ylabel('Difference [% of CTL]');
if ischar(return_period),if strcmpi(return_period,'AED'),ylim([-50 150]);end;end
hold on
errorbar(1:n_countries+1,res.climate_signal_mean,res.climate_signal_min,res.climate_signal_max,'.k');
for country_i=1:n_countries+1
    text(country_i-0.4,res.climate_signal_mean(country_i)+sign(res.climate_signal_mean(country_i))*5,sprintf('%i%%',ceil(res.climate_signal_mean(country_i))));
end % country_i

end % schwierz_etal_2010