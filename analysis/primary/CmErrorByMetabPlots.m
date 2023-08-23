% ------------------------------------------------------------
% Script:  CmErrorByMetabPlots.m
% Summary: renders a shaded curve indicating the range of possible
%          propagated error values for a particular metabolite of interest,
%          using the set of highest and lowest reported parameter
%          uncertainties from literature.
% 
% Usage:   >> CmErrorByMetabPlots
% 
% Input:   None. Script constants and uncertainties for this analysis are 
%          set in the following JSON file: ./data/CmErrorMetabs.json
%
% Output:  @file - PNG with filename: error-by-metabolite-[metabolite].png
%           that contains a plot of coefficient of variation (CV) across
%           reference error dAm (CRLB), as a shaded region bounded by the CV
%           derived from the highest and lowest reported set of parameter
%           uncertainties, as specified in CmErrorMetabs.json
%          @file - .MAT file containing metadata from the analysis,
%           including notes, results, and intermediary analysis variables. 
%
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
%
% Columbia University, 2023
% ------------------------------------------------------------
function CmErrorByMetabPlots(config)
    close all; clc;
    disp('Running Error Analysis by Metabolite...');

    % script parameters
    n        = config.CmErrorByMetabPlots.n;                % vector size
    notes    = config.CmErrorByMetabPlots.notes;            % analysis summary
    ds       = CmUtils.load_dataset('CmErrorMetabs.json');  % load relevant dataset
    meta     = config.CmErrorByMetabPlots.metabolites;      % metabolites to analyze
    meta_str = config.CmErrorByMetabPlots.meta_str;         % figure labels
    baseline = config.CmErrorByMetabPlots.baseline;         % error of reference (CRLB)
    refs     = config.CmErrorByMetabPlots.parameters;       % quantification parameters to include

    % display limits for metabolite as [xmax, ymax]
    axis_limits = {
       [50,50],[50,50],[50,50], ... % Cho, Cr, NAA
       [50,50],[70,70],[50,50]  ... % GABA, Gln, Glu
    };

    % parameter colors for shaded error ranges as HEX
    %   {blue, light blue, red, green, pink, orange}
    colors = {'#0076BA','#76D6FF','#EE220C','#1DB100','#FF42A1','#FEAE00'}';
    colors_rgb = {
        [0 118 186], [118 214 255], ... % blue, light blue
        [238 34 12], [29 177 0],    ... % red, green
        [255 66 161], [254 174 0]       % pink, orange
    };

    % set analysis options
    params_upper.options.metabolites = meta;   % calculate for all metabolites
    params_lower.options.metabolites = meta;
    params_upper.options.includedTerms = refs; % include all parameters in error
    params_lower.options.includedTerms = refs;

    % format dataset for error propagation
    params_upper.constants.params = ds.constants;                        % shared constants
    params_lower.constants.params = ds.constants;
    params_upper.constants.metabolites = ds.metabolites;                 % metab-specific constants
    params_lower.constants.metabolites = ds.metabolites;
    params_upper.errors.params = ds.errors.upper_bound.params;           % parameter errors
    params_lower.errors.params = ds.errors.lower_bound.params;
    params_upper.errors.metabolites = ds.errors.upper_bound.metabolites; % metabolite errors
    params_lower.errors.metabolites = ds.errors.lower_bound.metabolites;

    % covariances
    params_upper.covariances = ds.covariances.upper_bound;
    params_lower.covariances = ds.covariances.lower_bound;

    % sweep across dAm values
    for i = 1:length(meta)
        upper_dAm_temp = params_upper.constants.metabolites.(meta{i}).Am;
        lower_dAm_temp = params_lower.constants.metabolites.(meta{i}).Am;
        params_upper.errors.metabolites.(meta{i}).dAm = linspace(0,upper_dAm_temp,n);
        params_lower.errors.metabolites.(meta{i}).dAm = linspace(0,lower_dAm_temp,n);
    end

    % perform error propagation analysis
    [Cm_upper, dCm_upper, md_upper] = CmErrorByMetabolite(params_upper);
    [Cm_lower, dCm_lower, md_lower] = CmErrorByMetabolite(params_lower);

    % save raw results & metadata
    params.n = n; params.notes = notes;
    save(strcat(config.paths.res_dir,tempname('.'), ...
        '-error-by-metabolite','.mat'),             ...
        'md_upper','md_lower','params');

    % generate plots for each metabolite
    for i = 1:length(meta)

        % select metabolite
        field = meta{i};

        % relative error for each reference (dAm/Am*100)
        error_range_upper = params_upper.errors.metabolites.(field).dAm;
        error_range_lower = params_lower.errors.metabolites.(field).dAm;
        dref_upper = error_range_upper/params_upper.constants.metabolites.(field).Am*100;
        dref_lower = error_range_lower/params_lower.constants.metabolites.(field).Am*100;
        dCmPerc_upper = dCm_upper.(field)/Cm_upper.(field) * 100; 
        dCmPerc_lower = dCm_lower.(field)/Cm_lower.(field) * 100; 

        % initialize figure, set display parameters
        fig = figure; hold on; grid on;
        fig.Position = [100,100,500,350];
        set(fig, 'Color', [1 1 1]);

        % label title and axes
        title( meta_str{i},           ...
               'FontWeight','normal', ...
               'Fontsize',24 );
        ax = gca; 
        ax.FontSize = 20;
        ax.FontName = 'Calibri';
        xl = xlabel(strcat('\deltaS_m [%]'));
        yl = ylabel('CV, \deltaC_m/C_m [%]');
        set(xl,'FontSize',20);
        set(yl,'FontSize',20);
        xlim([0 axis_limits{i}(1)]);
        ylim([0 axis_limits{i}(2)]);

        % fill region between charts
        patch( [dref_upper fliplr(dref_lower)],       ...
               [dCmPerc_upper fliplr(dCmPerc_lower)], ...
               colors_rgb{i}./256,                    ...
               'FaceAlpha',.3 );
        plot( dref_upper,dCmPerc_upper, ...
              'color',colors{i},        ...
              'LineStyle','-',          ...
              'LineWidth',2 );
        plot( dref_lower,dCmPerc_lower, ...
              'color',colors{i},        ...
              'LineStyle','-',          ...
              'LineWidth',2 );

        % annotate with vertical line @ ref CRLB
        ref_crlb = ds.errors.upper_bound.metabolites.(field).dAm / ds.metabolites.(field).Am;
        xl = xline(ref_crlb*100,'--');
        xl.FontSize  = 18; 
        xl.LineWidth = 2;
        xl.FontName  = 'Calibri';
        xl.LabelHorizontalAlignment = 'left';
        xl.LabelVerticalAlignment   = 'top';

        % add linear term for reference
        lterm = linspace(0,100,n);
        plot( lterm,lterm,      ...
              'LineStyle','--', ...
              'color','black' );

        % annotate with error range
        dim = [.56 .21 .34 .17];
        if strcmp(field,'gaba')
            dim(1) = .61;
        end
        t1low  = strrep(compose("%.3f",params_lower.errors.metabolites.(field).dT1m),'0.','.');
        t1high = strrep(compose("%.3f",params_upper.errors.metabolites.(field).dT1m),'0.','.');
        t2low  = strrep(compose("%.3f",params_lower.errors.metabolites.(field).dT2m),'0.','.');
        t2high = strrep(compose("%.3f",params_upper.errors.metabolites.(field).dT2m),'0.','.');
        str = {
            strcat(t1low," s \leq ",'\deltaT1_m', " \leq ",t1high," s"), ...
            strcat(t2low," s \leq ",'\deltaT2_m', " \leq ",t2high," s")  ...
        };
        a = annotation(                    ...
            'textbox',dim,                 ...
            'String',str,                  ...
            'verticalalignment', 'middle', ...
            'horizontalalignment','center' );
        a.FontSize  = 15;
        a.LineStyle = 'None';
        a.BackgroundColor = '#F5F3F6';

        % save image
        saveas(fig,strcat(config.paths.res_dir,'error-by-metabolite-',field,'.png'));
    end
end