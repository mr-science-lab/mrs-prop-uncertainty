% ------------------------------------------------------------
% Script:  CmErrorByParamPlots.m
% Summary: performs an analysis of error on absolute quantification of MR
%          spectroscopy datasets by parameter. renders a shaded curve
%          indicated the range of possible contributions to the overall
%          propagated error by a single quantification parameter.
% 
% Usage:   >> CmErrorByParamPlots
%
% Input:   None. Script constants and uncertainties for this analysis are 
%          set in the following JSON file: ./data/CmErrorParams.json
%
% Output:  @file - PNG images with filename: single-ref-[param].png
%          @file - .MAT file containing analysis metadata, including including 
%           notes, results, and intermediary analysis variables. 
%
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
% 
% Columbia University, 2023
% ------------------------------------------------------------
function CmErrorByParamPlots(config)
    close all; clc; 

    % script parameters
    n     = config.CmErrorByParamPlots.n;               % vector size
    notes = config.CmErrorByParamPlots.notes;           % notes are included in output saved object
    ds    = CmUtils.load_dataset('CmErrorParams.json'); % quantification constants & uncertainties
    parameters = config.CmErrorByParamPlots.parameters; % single reference quantification parameters
    units = config.CmErrorByParamPlots.units;           % parameter to unit mapping

    % display limits for parameter as [xmax, ymax]
    axis_limits = {
        [5,10] , [10,20], [10,20],   ...  % Sw, T1m, T2m
        [10,20], [10,20],            ...  % f_i
        [1,2],   [1,2],              ...  % TE,TR
        [5,10],  [5,10],  [5,10]     ...  % T1_i
        [5,10],  [5,10],  [5,10]     ...  % T2_i
        [10,20], [10,20], [10,20],   ...  % cw_i
    };

    % parameter colors for shaded error ranges as HEX
    %   {grey, blue, light blue, green (x3), orange, purple,
    %   blue (x3), light blue (x3), red (x3)}
    colors = { '#929292','#0076BA','#76D6FF','#1DB100','#1DB100',           ...
               '#FEAE00','#9830A0','#0076BA','#0076BA','#0076BA','#76D6FF', ...
               '#76D6FF','#76D6FF','#EE220C','#EE220C','#EE220C'}';
    colors_rgb = {
        [146 146 146], [0 118 186],   [118 214 255], ...  % grey, blue, light blue
        [29 177 0],    [29 177 0],                   ...  % green, green
        [254 174 0],   [152 48 160],                 ...  % pink, orange, purple
        [0 118 186],   [0 118 186],   [0 118 186],   ...  % blue (T1_i)
        [118 214 255], [118 214 255], [118 214 255], ...  % light blue (T2_i)
        [238 34 12],   [238 34 12],   [238 34 12]    ...  % red (cw_i)
    };

    % format dataset for error propagation
    upper.constants = ds.constants;                                     % upper bound shared constants
    lower.constants = ds.constants;                                     % lower bound shared constants
    upper.errors    = ds.errors.upper_bound;                            % upper bound uncertainties
    lower.errors    = ds.errors.lower_bound;                            % lower bound uncertainties
    upper.options.reference_var = 'Sm';                                 % set linear term (reference)
    lower.options.reference_var = 'Sm';
    upper.options.error_range   = linspace(0,upper.constants.Sm, n);    % calculate over ref uncertainty
    lower.options.error_range   = linspace(0,lower.constants.Sm, n);
    upper.covariances = ds.covariances.upper_bound;                     % upper bound covariances
    lower.covariances = ds.covariances.lower_bound;                     % lower bound covariances

    % generate a figure for each parameter
    for i = 1:length(parameters)

        % set current parameter for analysis
        param = parameters{i};
        upper.options.includedTerms = {'Sm', param}; % linear term + variable of interest 
        lower.options.includedTerms = {'Sm', param}; 

        % calculate propagated error (dCm)
        disp(strcat('Running Uncertainty Analysis on Parameter:',param));
        [Cm_upper, dCm_upper, md_upper] = CmErrorByParameter(upper);
        [Cm_lower, dCm_lower, md_lower] = CmErrorByParameter(lower);

        % relative error of each reference (dSm/Sm*100), same for both
        dref_upper = upper.options.error_range / upper.constants.Sm*100;
        dref_lower = lower.options.error_range / lower.constants.Sm*100;

        % normalized coefficient of variation, estimate
        dCmrel_upper = dCm_upper / Cm_upper * 100;
        dCmrel_lower = dCm_lower / Cm_lower * 100;

        % format and style figure
        fig = figure; hold on; grid on; 
        fig.Position = [100,100,500,350];                   % set size, position
        set(fig,'Color',[1 1 1]);                           % set figure background
        param_str = CmUtils.format_param_as_string(param);  % format string for title
        title(strcat({'\delta'},param_str), ...             % add & style figure title
            'FontWeight','normal',          ...
            'Fontsize',14,                  ...
            'FontName','Calibri');
        xl = xlabel('\deltaS_m [%]');                       % x-axis label
        yl = ylabel('CV, \deltaC_m/C_m [%]');               % y-axis label

        % set font sizes, axis limits
        ax = gca;
        ax.FontSize = 20;
        ax.FontName = 'Calibri';
        set(xl,'FontSize',20);
        set(yl,'FontSize',22);
        set(xl,'FontName','Calibri');
        set(yl,'FontName','Calibri');
        set(yl,'color',colors{i});
        xlim([0,axis_limits{i}(1)]);
        ylim([0,axis_limits{i}(2)]);

        % fill region between charts
        patch([dref_upper fliplr(dref_lower)],       ...
              [dCmrel_upper fliplr(dCmrel_lower)],   ...
              colors_rgb{i}./256,'FaceAlpha',.3);

        % render boundaries
        plot(dref_upper,dCmrel_upper,'color',colors{i},'LineStyle','-','LineWidth',2);
        plot(dref_lower,dCmrel_lower,'color',colors{i},'LineStyle','-','LineWidth',2);

        % add annotation with error range
        dim = [.46 .76 .41 .1];
        if strcmp('TE',param) || strcmp('TR',param)
            formatspec = "%.e";
        else
            formatspec = "%.3f";
        end
        low  = strrep(compose(formatspec,lower.errors.(strcat('d',param))),'0.','.');
        high = strrep(compose(formatspec,upper.errors.(strcat('d',param))),'0.','.');
        str = strcat(low,{' '},units.(param)," \leq ",'\delta', ...
            CmUtils.format_param_as_string(param)," \leq ",high,{' '},units.(param));
        a = annotation('textbox',dim,'String',str,'verticalalignment', ...
            'middle','horizontalalignment','center');
        a.FontSize = 16;
        a.FontName = 'Calibri';
        a.LineStyle = 'None';
        a.BackgroundColor = '#F5F3F6';
        if strcmp(param,'cw_white') % decrease font size for cw_white
            a.FontSize = 15;        % ensure legend fits in textbox
        end

        % generate linear term for reference
        % (i.e. propagated error as a result of a single linear term, dSm)
        reference.constants = ds.constants;          % use same constants as previous analysis
        reference.errors    = ds.errors.upper_bound; % placeholder; doesn't matter since terms are not included
        reference.options.reference_var = 'Sm';      % error of reference (CRLB)
        reference.options.error_range   = linspace(0,upper.constants.Sm, n);
        reference.options.includedTerms = {'Sm'};    % do not include any terms other than CRLB
        reference.covariances = ds.covariances.upper_bound;
        [Cm_ref, dCm_ref] = CmErrorByParameter(reference);
        dreference = reference.options.error_range / reference.constants.Sm*100;
        dCmrel_reference = dCm_ref / Cm_ref * 100;

        % include linear term to figure as dashed line
        plot( dreference,dCmrel_reference, ...
              'color',colors{i},           ...
              'LineStyle','--',            ...
              'LineWidth',2 );

        % determine percent contribution of parameter compared to CRLB
        % -- ratio of parameter to overall uncertainty
        % -- (total unceratinty - CRLB uncertainty) / total uncertainty
        mean_unceratinty  = (dCmrel_upper+dCmrel_lower)/2;
        param_percent_avg = (mean_unceratinty - dCmrel_reference) ./ mean_unceratinty;

        % add percent influence curve to figure
        yyaxis right;
        plot( dreference, param_percent_avg, ...
              'color', 'black',              ...
              'LineWidth', 2);
        ax.YAxis(2).Color = 'k';
        ax.YAxis(1).Color = colors{i};
        yl = ylabel('Parameter Influence');
        set(yl,'FontSize',20);

        % save image as PNG
        exportgraphics(fig,strcat(config.paths.res_dir,'single-ref-',param,'.png'),'Resolution',1000);

        % save analysis parameters as MAT
        params.n = n; 
        params.targetparam = param; 
        params.notes = notes;
        save(strcat(config.paths.res_dir,tempname('.'),'-',param,'.mat'), ...
            'md_upper','md_lower','params');
    end
end