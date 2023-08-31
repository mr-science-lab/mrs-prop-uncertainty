% ------------------------------------------------------------
% Script:  CmErrorParamInfluence.m
% Summary: Plots partial derviative terms over T1m, T2m, TE, TR and fi_c.
%          Examines the effect of sequence and metabolite-specific terms in
%          the analytical solution to partial derivatives on the overall
%          propagated uncertainty. 
%
% Output:  @file - two PNG images containing the partial derviatives as a 
%          T1m and T2m values, with annotated vertical lines indicating
%          reported T1, T2 values for various metabolites.
%
% Usage:   >> CmErrorParamInfluence
%
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
%
% Columbia University, 2023
% ------------------------------------------------------------
function CmErrorParamInfluence(config)
    close all; clc;
    disp('Plotting analytical solution to partial derviatives by parameter...');

    % generate single plot with all metabolites
    % separate individual plots with shaded region
    plots = {
        {'cho','cr','gaba','gln','glu','naa'}, ...
        {'cho'}, {'cr'}, {'gaba'}, {'gln'}, {'glu'}, {'naa'}
    };
    colors.cho  = '#0076BA';    colors_rgb.cho  = [0 118 186];
    colors.cr   = '#76D6FF';    colors_rgb.cr   = [118 214 255];
    colors.gaba = '#1DB100';    colors_rgb.gaba = [29 177 0];
    colors.gln  = '#FF42A1';    colors_rgb.gln  = [255 66 161];
    colors.glu  = '#FEAE00';    colors_rgb.glu  = [254 174 0];
    colors.naa  = '#EE220C';    colors_rgb.naa  = [238 34 12];

    for ploti = 1:length(plots)

        % TERM: T1m
        n  = config.CmErrorParamInfluence.n;  % vector size
        ds = CmUtils.load_dataset('CmErrorParams.json');

        % partial of T1m
        T1m  = linspace(0.1,1.7,n);
        dT1m = zeros(1,n);
        for i = 1:n
            ds.constants.T1m = T1m(i);
            CmObj = CmPartials(ds.constants);
            dT1m(i) = CmObj.partial_T1m;
        end

        % plot T1m curve
        fig = figure; 
        fig.Position = [1000,1000,850,230];
        if ploti ~= 1
            fig.Position = [1000,1000,425,230];
        end
        plot( T1m,dT1m/dT1m(1), ...
              'LineWidth',2,    ...
              'Color','black');
        ax = gca;
        grid on;
        if ploti == 1
            xl = xlabel('{\it T}_{1m} (s)');
            yl = ylabel('\partial{C_m}/\partial{{\it T}_{1m}}', ...
                'FontWeight','normal');
        end
        xticks(linspace(0.1,1.7,17));
        xlim([0.1,1.7]);
        if ploti ~= 1
            xlim([0.7,1.7]);
        end
        set(gca,'GridLineStyle','-');
        set(gca,'yticklabel',[]);     % exclude y-axis tick labels
        set(gca,'ytick',[]);          % exclude y-axis tick marks
        set(gca,'XMinorTick','on');   % include minor x-axis tick marks
        set(gca,'xminorgrid','on');   % include minor grid lines
        set(fig,'Color',[1 1 1]);     % set white figure background
        ax.FontSize = 18;             % axes label font size
        ax.FontName = 'Calibri';      % axes label font name

        % annotate with T1 values of metabolites, vertical lines. 
        ms = CmUtils.load_dataset('CmErrorMetabs.json');

        metas = plots{ploti};
        for i = 1:length(metas)

            % mark T1 of metabolite as vertical line
            linestyle = '-';
            color = 'black';
            if ploti == 1
                color = colors.(metas{i});
            end
            xl = xline( ms.metabolites.(metas{i}).T1m,                        ...
                        linestyle,{CmUtils.format_metab_as_string(metas{i})}, ...
                        'color',color );
            xl.FontSize  = 18;
            xl.LineWidth = 2;
            xl.FontName  = 'Calibri';
            xl.LabelHorizontalAlignment = 'left';
            if strcmp(metas{i},'gaba') || strcmp(metas{i},'gln') || strcmp(metas{i},'cr')
                xl.LabelHorizontalAlignment = 'right';
            end
            xl.LabelVerticalAlignment = 'bottom';
            if strcmp(metas{i},'gln') && ploti == 1
                xl.LabelHorizontalAlignment = ' center';
                xl.FontSize = 16;
            end
            if strcmp(metas{i},'glu') && ploti == 1
                xl.LabelHorizontalAlignment = 'center';
                xl.FontSize = 16;
            end

            if ploti ~= 1
               % use to mark borders of shaded region
               ylimit = max(dT1m/dT1m(1));

               % add shaded region for reported error range, use lower bound
               interval = ms.errors.lower_bound.metabolites.(metas{i}).dT1m;
               lower    = ms.metabolites.(metas{i}).T1m-interval; 
               upper    = ms.metabolites.(metas{i}).T1m+interval;
               ylimit   = max(dT1m/dT1m(1));
               patch([lower,lower,upper,upper],[0 ylimit ylimit 0], ...
                     colors_rgb.(metas{i})./256,                    ...
                     'FaceAlpha',.3,                                ...
                     'LineStyle','none');

               % add shaded region for upper bound (lower opacity)
               interval = ms.errors.upper_bound.metabolites.(metas{i}).dT1m;
               lower    = ms.metabolites.(metas{i}).T1m-interval; 
               upper    = ms.metabolites.(metas{i}).T1m+interval;
               patch([lower,lower,upper,upper],[0 ylimit ylimit 0], ...
                     colors_rgb.(metas{i})./256,                    ...
                     'FaceAlpha',.15,                               ...
                     'LineStyle','none');
            end
        end

        % save image
        exportgraphics(fig,strcat(config.paths.res_dir,'parameter-partial-T1-',int2str(ploti),'.png'),'Resolution',2000);

        % TERM: T2m
        clear c T1m;
        ds = CmUtils.load_dataset('CmErrorParams.json');

        % partial of T2m
        T2m  = linspace(0.01,1,n);
        dT2m = zeros(1,n);
        for i = 1:n
            ds.constants.T2m = T2m(i);
            CmObj = CmPartials(ds.constants);
            dT2m(i) = CmObj.partial_T2m;
        end

        % plot T2m curve
        fig = figure; 
        fig.Position = [1000,1000,850,230];
        if ploti ~= 1
            fig.Position = [1000,1000,425,230];
        end
        plot(T2m,abs(dT2m)/abs(dT2m(1)), ...
             'LineWidth',2,              ...
             'Color','black');
        ax = gca;
        grid on; 
        if ploti == 1
            xl = xlabel('{\it T}_{2m} (s)');
            yl = ylabel('|\partial{C_m}/\partial{{\it T}_{2m}}|', ...
                        'FontWeight','normal');
        end
        xticks(linspace(0.05,1,20));
        xlim([0.05,0.4]);
        ylim([0,.004]);
        ax.FontSize = 18;
        ax.FontName = 'Calibri';
        set(gca,'yticklabel',[]);         % exclude y-axis tick labels
        set(gca,'ytick',[]);              % exclude y-axis tick marks
        set(gca,'XMinorTick','on');       % include minor x-axis tick marks
        set(gca,'xminorgrid','on');       % include minor grid lines   
        set(fig,'Color',[1 1 1]);         % set white figure background

        % add error bars as shaded regions. 
        for i = 1:length(metas)
            color = 'black';
            if ploti == 1
                color = colors.(metas{i});
            end

            % mark T2 of metabolite as vertical line
            xl = xline(ms.metabolites.(metas{i}).T2m,'-', ...
                {CmUtils.format_metab_as_string(metas{i})},'color',color);
            xl.FontSize = 18; xl.LineWidth = 2; xl.FontName = 'Calibri';
            xl.LabelHorizontalAlignment = 'left';
            xl.LabelVerticalAlignment = 'top';
            if (strcmp(metas{i},'gln') || strcmp(metas{i},'gaba') || strcmp(metas{i},'glu'))
                if ploti == 1
                    xl.LabelHorizontalAlignment = 'center';
                else
                    xl.LabelHorizontalAlignment = 'right';
                end
            end

            if ploti ~= 1

               % add shaded region for reported error range, use lower bound
               interval = ms.errors.lower_bound.metabolites.(metas{i}).dT2m;
               lower    = ms.metabolites.(metas{i}).T2m-interval; 
               upper    = ms.metabolites.(metas{i}).T2m+interval;
               patch([lower,lower,upper,upper],[0 1 1 0], ...
                     colors_rgb.(metas{i})./256,          ...
                     'FaceAlpha',.3,                      ...
                     'LineStyle','none');

               % add shaded region for upper bound (lower opacity)
               interval = ms.errors.upper_bound.metabolites.(metas{i}).dT2m;
               lower    = ms.metabolites.(metas{i}).T2m-interval; 
               upper    = ms.metabolites.(metas{i}).T2m+interval;
               patch([lower,lower,upper,upper],[0 ylimit ylimit 0], ...
                     colors_rgb.(metas{i})./256,                    ...
                     'FaceAlpha',.15,                               ...
                     'LineStyle','none');
            end
        end

        % save image
        exportgraphics(fig,strcat(config.paths.res_dir,'parameter-partial-T2-',int2str(ploti),'.png'),'Resolution',2000);
    end
end