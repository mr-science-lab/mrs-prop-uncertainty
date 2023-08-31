% ------------------------------------------------------------
% Script:  CmErrorSensitivityPlots
% Summary: Plots the effect of f_i covariances on propagated uncertainty.
%          A range of covariances from 0 to std(X)*std(Y) is plotted, where
%          X,Y are f_grey and f_white. The lower bound uncertainty of f_i
%          is chosen for this analysis.
%
% Usage:   >> CmErrorSensitivityPlots
%
% Input:   None.
%
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
% 
% Columbia University, 2023
% ------------------------------------------------------------
function CmErrorSensitivityPlots(config)
    close all; clc; 

    % script parameters
    n  = config.CmErrorSensitivityPlots.n;                 % vector size (uncertainty curves)
    ns = config.CmErrorSensitivityPlots.ns;                % vector size (sensitivity curves)
    ds = CmUtils.load_dataset('CmErrorMetabs.json');       % metabolite parameters and constants

    % reference parameters
    meta     = config.CmErrorSensitivityPlots.metabolites; % metabolites
    meta_str = config.CmErrorSensitivityPlots.meta_str;    % metabolite string labels for figures
    refs     = config.CmErrorSensitivityPlots.parameters;  % include all parameters for tissue-specific internal water reference

    % set display parameters
    % axes limits [xmin,xmax,ymin,ymax]
    axis_limits = {[0,15,0,15],[0,15,0,15],[0,15,0,15],[0,15,0,15],[0,15,0,15],[0,15,0,15]};

    % individual metabolite colors
    % {blue, light blue, red, green, pink, orange}
    colors = {'#0076BA','#76D6FF','#EE220C',...
        '#1DB100','#FF42A1','#FEAE00'}';        % HEX color values
    colors_rgb = {                              % include RGB of each color
        [0 118 186], [118 214 255], ...         % blue, light blue
        [238 34 12], [29 177 0],    ...         % red, green
        [255 66 161], [254 174 0]               % pink, orange
    };

    % set analysis options
    params_upper.options.metabolites   = meta;                  % calculate for all metabolites
    params_lower.options.metabolites   = meta;
    params_upper.options.includedTerms = refs;                  % include all parameters in error
    params_lower.options.includedTerms = refs;

    % format dataset for error propagation
    params_upper.constants.params = ds.constants;               % shared constants
    params_lower.constants.params = ds.constants;
    params_upper.constants.metabolites = ds.metabolites;        % metab-specific constants
    params_lower.constants.metabolites = ds.metabolites;
    params_upper.errors.params = ds.errors.upper_bound.params;  % parameter errors
    params_lower.errors.params = ds.errors.lower_bound.params;
    params_upper.errors.metabolites = ds.errors.upper_bound.metabolites;
    params_lower.errors.metabolites = ds.errors.lower_bound.metabolites;

    % set range of covariances
    f_covs = linspace(-0.02^2, 0.02^2, ns);

    % run for each metabolite
    for meta_index = 1:length(meta)

        % track propagated uncertainty at CRLB 0
        uncert_at_crlb0{meta_index} = [];

        % calculate for single metabolite
        params_lower.options.metabolites = {meta{meta_index}};

        % generate a single total unceratinty for each covariance
        for j = 1:length(f_covs)

            % covariances
            params_lower.covariances = ds.covariances.lower_bound;

            % set covariance based on iteration
            disp(['covariance ' num2str(j) ' of ' num2str(length(f_covs))]);
            disp(['cov(f_white,f_grey):' num2str(f_covs(j))]);
            params_lower.covariances.f_grey_f_white = f_covs(j);

            % sweep across dSm values
            lower_dSm_temp = params_lower.constants.metabolites.(meta{meta_index}).Sm;
            params_lower.errors.metabolites.(meta{meta_index}).dSm = linspace(0,lower_dSm_temp,n);

            % perform error propagation analysis
            [Cm_lower, dCm_lower, md_lower] = CmErrorByMetabolite(params_lower);

            % select metabolite
            field = meta{meta_index};

            % relative error for each reference (dSm/Sm*100)
            error_range_lower = params_lower.errors.metabolites.(field).dSm;
            dref_lower = error_range_lower/params_lower.constants.metabolites.(field).Sm*100;
            dCmPerc_lower = dCm_lower.(field)/Cm_lower.(field) * 100; 

            % initialize figure, set display parameters
            if (j == 1)
                fig = figure; hold on; grid on;
                fig.Position = [100,100,500,350];
                set(fig, 'Color', [1 1 1]);
            end

            % label title and axes
            title(meta_str{meta_index},'FontWeight','normal','Fontsize',24);
            ax = gca;
            ax.FontSize = 20;
            ax.FontName = 'Calibri';
            xl = xlabel(strcat('\deltaS_m [%]'));
            yl = ylabel('CV, \deltaC_m/C_m [%]');
            set(xl,'FontSize',20);
            set(yl,'FontSize',20);
            xlim([axis_limits{meta_index}(1) axis_limits{meta_index}(2)]);
            ylim([axis_limits{meta_index}(3) axis_limits{1}(4)]);

            % include 3 distinct curves
            if (f_covs(j) == 0)
                % no covariance (black)
                plot( dref_lower,dCmPerc_lower, ...
                      'color','black',          ...
                      'LineStyle','-',          ...
                      'LineWidth',1.5 );
            elseif (j == 1)
                % lowest cov
                plot( dref_lower,dCmPerc_lower,   ...
                      'color',colors{meta_index}, ...
                      'LineStyle','-',            ...
                      'LineWidth',1.5 );
            elseif (j == length(f_covs))
                % highest cov
                plot( dref_lower,dCmPerc_lower,   ...
                      'color',colors{meta_index}, ...
                      'LineStyle','-',            ...
                      'LineWidth',1.5);
            end

            % save the propagated uncertainties @ CRLB = 0
            uncert_at_crlb0{meta_index} = [uncert_at_crlb0{meta_index} dCmPerc_lower(1)];

        end    

        % add linear term for reference
        lterm = linspace(0,100,n);
        plot(lterm,lterm,'LineStyle','--','color','black','LineWidth',1);

        % save image
        exportgraphics(fig,strcat(config.paths.res_dir,'sensitivity-by-metabolite-',field,'.png'),'Resolution',2000);
        saveas(fig,strcat(config.paths.res_dir,'sensitivity-by-metabolite-',field,'.fig'));

    end

    % plot sensitivity as a function of covariance
    fig = figure; hold on; grid on;             % all metabolites on single figure
    fig.Position = [100,100,500,350];           % figure dimensions
    set(fig, 'Color', [1 1 1]);                 % set background color to white
    for plot_index = 1:6
        plot( f_covs,uncert_at_crlb0{plot_index}, ...
              'color',colors{plot_index},         ...
              'LineStyle','-',                    ...
              'LineWidth',2 );
        hold on; 
    end
    legend('Cho','Cr','NAA','GABA','Gln','Glu','Location','southeast');
    xlim([-0.02^2, 0.02^2]);                    % x-axis limits
    xl = xlabel('cov(f_{grey},f_{white})');     % y-axis label
    yl = ylabel('CV, \deltaC_m/C_m [%]');       % x-axis label
    set(xl,'FontSize',20);                      % set x-axis font size
    set(yl,'FontSize',20);                      % set y-axis font size
    ax = gca; 
    ax.FontSize = 16; 
    ax.FontName = 'Calibri';
    
    % save image
    exportgraphics(fig,strcat(config.paths.res_dir,'sensitivity-curves','.png'),'Resolution',2000);
    saveas(fig,strcat(config.paths.res_dir,'sensitivity-curves','.fig'));
end
