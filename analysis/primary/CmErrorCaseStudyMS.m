% ------------------------------------------------------------
% Script:  CmErrorCaseStudyMS.m
% Summary: Runs a Monte Carlo analysis to estimate the propagated
%          uncertainty of absolute quantification of subjects in 
%          the following MS study
%
% Inputs:  Quantification parameters are organized by subject, and
%          can be found here: 
%          ./data/ms_config/*
%
% Usage:   >> CmErrorCaseStudyMS.m
%
% Reference: 
%          Swanberg KM, Prinsen H, DeStefano K, et al. In vivo evidence 
%          of differential frontal cortex metabolic abnormalities in progressive
%          and relapsing-remitting multiple sclerosis. NMR Biomed. 2021;34(11).
%          doi:10.1002/nbm.4590
%
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
%
% Columbia University, 2023
% ------------------------------------------------------------
function CmErrorCaseStudyMS(config)
    close all; clc;
    disp('Running Case Study analysis: MS dataset...');

    % analysis parameters
	n     = config.CmErrorCaseStudyMS.n;     % number of simulated Cm Values
    nbins = config.CmErrorCaseStudyMS.nbins; % number of bins in simulated histogram
    ds    = CmUtils.load_dataset(config.CmErrorCaseStudyMS.filename);

    %% Perform Monte Carlo simulation on metabolite: NAA
    disp('simulating metabolite: NAA');

    % include all parameters in propagated uncertainty estimate
    all_refs = { ...
        { 'Am' }, ...
        { 'Am', 'S_tCr' }, ...
        { 'Am', 'S_tCr', 'T1_tCr' }, ...
        { 'Am', 'S_tCr', 'T1_tCr', 'T2_tCr' }, ...
        { 'Am', 'S_tCr', 'T1_tCr', 'T2_tCr', 'T1m' }, ...
        { 'Am', 'S_tCr', 'T1_tCr', 'T2_tCr', 'T1m', 'T2m' }, ...
        { 'Am', 'S_tCr', 'T1_tCr', 'T2_tCr', 'T1m', 'T2m', 'TE', 'TR' } ...
    };

    % initialize figure
    colors = {'#5E5E5E','#5E5E5E','#5E5E5E','#5E5E5E','#5E5E5E','#5E5E5E','#0076BA'};
    fig = figure;                     % initialize figure
    fig.Position = [100,100,500,400]; % figure dimensions, size
    set(fig,'Color',[1 1 1]);         % set white figure background

    % simulate each subset of parameters
    for refset_index = 1:length(all_refs)

        % display current parameter of interest
        disp(strcat('Parameter of interest: ',all_refs{refset_index}{end}));

        % load quantification parameters
        quant_params = {};
        quant_params.Am      = ds.constants.metabolites.tNAA.Am;
        quant_params.T1m     = ds.constants.metabolites.tNAA.T1;
        quant_params.T2m     = ds.constants.metabolites.tNAA.T2;
        quant_params.S_tCr   = ds.constants.metabolites.tCr.Am;
        quant_params.T1_tCr  = ds.constants.metabolites.tCr.T1;  % use T1 of Cr for tCr
        quant_params.T2_tCr  = ds.constants.metabolites.tCr.T2;  % use T2 of Cr for tCr
        quant_params.TE      = ds.constants.acquisition.TE;
        quant_params.TR      = ds.constants.acquisition.TR;
        quant_params.C_tCr   = ds.ref_concentration;             % 10 mM internal reference

        % set quantification uncertainties
        quant_uncert.dAm      = ds.uncertainties.metabolites.tNAA.dAm;
        quant_uncert.dT1m     = ds.uncertainties.metabolites.tNAA.dT1;
        quant_uncert.dT2m     = ds.uncertainties.metabolites.tNAA.dT2;
        quant_uncert.dS_tCr   = ds.uncertainties.metabolites.tCr.dAm;
        quant_uncert.dT1_tCr  = ds.uncertainties.metabolites.tCr.dT1;
        quant_uncert.dT2_tCr  = ds.uncertainties.metabolites.tCr.dT2;
        quant_uncert.dTE      = ds.uncertainties.acquisition.dTE;
        quant_uncert.dTR      = ds.uncertainties.acquisition.dTR;

        % set Monte Carlo simulation parameters
        sim_params = quant_params;      % handle constants with no variance in simulation
        refs = all_refs{refset_index};  % set current set of simulated parameters

        % simulate a distribution for each reference parameter
        for ref_index = 1:length(refs)

            % set current parameter
            curr_ref = refs{ref_index}; 

            % set lower bound of truncation based on param
            if ismember(curr_ref,{'T1m','T2m','T1_tCr','T2_tCr','TE','TR'})
                truncate_lower = 0.001; % min time of 1 ms
            else
                truncate_lower = 0; % default to non-negative A.U. values
            end

            % simulate parameter distribution
            mu = quant_params.(curr_ref);
            sigma = quant_uncert.(strcat('d',curr_ref));
            pd = makedist('normal','mu',mu,'sigma',sigma);
            pd_truncated = truncate(pd,truncate_lower,inf); % truncate distribution
            pd_random = random(pd_truncated,1,n);           % generate random values

            % replace simulation parameter constant with simulated distribution
            sim_params.(curr_ref) = pd_random;

        end

        % simulate absolute quantification (in mM, since tCr is in mM)
        Cm = CmUtils.absolute_quantification_tCr(sim_params);
        if (refset_index == 1)
            Cm_crlb = Cm;
        end
        Cm(Cm > 50) = []; % outlier exclusion criteria

        % plot best fit line
        h = histfit(Cm,nbins,'kernel');
        h(1).Visible = 'off';
        h(2).Visible = 'off';
        x = h(2).XData; % density curve
        y = h(2).YData;
        y_orig = y;
        if (refset_index == 1)
            y_orig_crlb = y_orig;
        end
        y     = y / max(y); % normalize y-axis
        x_bar = h(1).XData;
        y_bar = h(1).YData / max(h(1).YData);
        lw    = 2;         % linewidth
        if refset_index == length(all_refs)
            lw = 4; 
        end
        plot(x,y,                          ...
             'Color',colors{refset_index}, ...
             'LineWidth',1.5);
        hold on; 
        grid on;

        % plot histogram only for CRLB distribution
        if refset_index == 1

            % generate normalized histogram data
            [counts, edges] = histcounts(Cm,nbins);
            counts = counts / max(y_orig);

            % styling parameters
            bar(edges(1:end-1),counts,  ...
                'FaceColor', '#0076BA', ...
                'FaceAlpha', 0.8,       ...
                'EdgeAlpha',0);
            axis square;

            % axes labels
            title('\rm tNAA',    ...
                  'FontSize',24, ...
                  'FontName','Calibri');
            xlabel('C_m (mM)'); 
            ylabel('Frequency');
            ax = gca; 
            ax.FontSize = 18;
            ax.FontName = 'Calibri';
            set(gca,'yticklabel',[])
            set(gca,'ytick',[])

          % set axis limit based on file
          ylim([0 1]);        % y-axis limits
          xlim([9.941 10.7]); % x-axis limits
                              % dataset: 01_highFWHM_M_HC  [9.941 10.7]
                              % dataset: 02_lowFWHM_F_RRMS [13.26 14.1]
        end

        % calculate CV_adjusted
        if refset_index == length(all_refs) || refset_index == 1
            CV_normal   = std(Cm)/mean(Cm)
            CV_adjusted = sqrt(exp(std(log(Cm))^2)-1)
        end

        % refresh figure, pause for user input
        drawnow;
    end

    % save the figure as PNG
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMS.filename,'_tNAA_param.png'));
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMS.filename,'_tNAA_param.fig'),'fig');

    % --- show full width distribution --- 
    % initialize figure
    fig = figure; 
    fig.Position = [100,100,500,400]; % position, dimensions
    set(fig,'Color',[1 1 1]);         % set white figure background

    % distribution with CRLB of Am only
    [counts, edges] = histcounts(Cm_crlb,nbins);
    counts = counts / max(y_orig_crlb);
    bar( edges(1:end-1),counts, ...
         'FaceColor','#0076BA', ...
         'FaceAlpha',0.9,       ...
         'EdgeAlpha',0);
    hold on; 

    % distribution of full propagated uncertainty
    [counts, edges] = histcounts(Cm,nbins);
    counts = counts / max(y_orig);
    bar( edges(1:end-1),counts, ...
         'FaceColor','#0076BA', ...
         'FaceAlpha',0.4,       ...
         'EdgeAlpha',0);

    % style figure
    title('\rm tNAA',    ...
          'FontSize',24, ...
          'FontName','Calibri'); 
    axis square;
    xlabel('C_m (mM)');
    ylabel('Frequency');
    ax = gca;
    ax.FontSize = 18;
    ax.FontName = 'Calibri';
    xlim([9 14.2]); % set x-axis limits based on dataset
    ylim([0 1]);    % normalize y-axis limits

    % determine best fit lognormal line
    h2 = histfit(Cm,nbins,'lognormal');
    h2(1).Visible = 'off';
    h2(2).Visible = 'off';
    x_lognormal   = h2(2).XData;                    % density curve
    y_lognormal   = h2(2).YData;
    y_lognormal   = y_lognormal / max(y_lognormal); % normalize y-axis

    % plot best fit line
    plot(x,y,'Color','#0076BA','LineWidth',2); hold on;                 % (kernel)
    plot(x_lognormal, y_lognormal, '--', 'Color', 'k', 'Linewidth', 1); % (lognormal)

    % save the figure as PNG
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMS.filename,'_tNAA_full.png'));
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMS.filename,'_tNAA_full.fig'),'fig');


    %% Perform Monte Carlo simulation on metabolite: Glu
    disp('simulating metabolite: Glu');

    % include all parameters in propagated uncertainty estimate
    all_refs = { ...
        { 'Am' }, ...
        { 'Am', 'S_tCr' }, ...
        { 'Am', 'S_tCr', 'T1_tCr' }, ...
        { 'Am', 'S_tCr', 'T1_tCr', 'T2_tCr' }, ...
        { 'Am', 'S_tCr', 'T1_tCr', 'T2_tCr', 'T1m' }, ...
        { 'Am', 'S_tCr', 'T1_tCr', 'T2_tCr', 'T1m', 'T2m' }, ...
        { 'Am', 'S_tCr', 'T1_tCr', 'T2_tCr', 'T1m', 'T2m', 'TE', 'TR' } ...
    };

    % set colors, figure parameters
    colors = {'#5E5E5E','#5E5E5E','#5E5E5E','#5E5E5E','#5E5E5E','#5E5E5E','#FEAE00'};
    fig = figure;                     % initialize figure
    fig.Position = [100,100,500,400]; % size, dimensions
    set(fig,'Color',[1 1 1]);         % set white figure background

    % simulate each subset of parameters
    for refset_index = 1:length(all_refs)

        % load quantification parameters
        quant_params = {};
        quant_params.Am      = ds.constants.metabolites.Glu.Am;
        quant_params.T1m     = ds.constants.metabolites.Glu.T1;
        quant_params.T2m     = ds.constants.metabolites.Glu.T2;
        quant_params.S_tCr   = ds.constants.metabolites.tCr.Am;
        quant_params.T1_tCr  = ds.constants.metabolites.tCr.T1;  % use T1 of Cr for tCr
        quant_params.T2_tCr  = ds.constants.metabolites.tCr.T2;  % use T2 of Cr for tCr
        quant_params.TE      = ds.constants.acquisition.TE;
        quant_params.TR      = ds.constants.acquisition.TR;
        quant_params.C_tCr   = ds.ref_concentration;             % 10 mM internal reference

        % set quantification uncertainties
        quant_uncert.dAm      = ds.uncertainties.metabolites.Glu.dAm;
        quant_uncert.dT1m     = ds.uncertainties.metabolites.Glu.dT1;
        quant_uncert.dT2m     = ds.uncertainties.metabolites.Glu.dT2;
        quant_uncert.dS_tCr   = ds.uncertainties.metabolites.tCr.dAm;
        quant_uncert.dT1_tCr  = ds.uncertainties.metabolites.tCr.dT1;
        quant_uncert.dT2_tCr  = ds.uncertainties.metabolites.tCr.dT2;
        quant_uncert.dTE      = ds.uncertainties.acquisition.dTE;
        quant_uncert.dTR      = ds.uncertainties.acquisition.dTR;

        % set Monte Carlo simulation parameters
        sim_params = quant_params;      % handle constants with no variance in simulation
        refs = all_refs{refset_index};  % set current set of simulated parameters

        % simulate a distribution for each reference parameter
        for ref_index = 1:length(refs)

            % set current parameter
            curr_ref = refs{ref_index}; 

            % set lower bound of truncation based on param
            if ismember(curr_ref,{'T1m','T2m','T1_tCr','T2_tCr','TE','TR'})
                truncate_lower = 0.001; % min time of 1 ms
            else
                truncate_lower = 0; % default to non-negative A.U. values
            end

            % simulate parameter distribution
            mu    = quant_params.(curr_ref);
            sigma = quant_uncert.(strcat('d',curr_ref));
            pd    = makedist('normal','mu',mu,'sigma',sigma);
            pd_truncated = truncate(pd,truncate_lower,inf); % set boundaries
            pd_random    = random(pd_truncated,1,n); % generate random values

            % replace simulation parameter constant with simulated distribution
            sim_params.(curr_ref) = pd_random;
        end

        % simulate absolute quantification (in mM, since tCr is in mM)
        Cm = CmUtils.absolute_quantification_tCr(sim_params);
        if (refset_index == 1)
            Cm_crlb = Cm;
        end
        Cm(Cm > 50) = []; % outlier exclusion criteria

        % plot best fit line
        h = histfit(Cm,nbins,'kernel');
        h(1).Visible = 'off';
        h(2).Visible = 'off';
        x = h(2).XData;
        y = h(2).YData;
        y_orig = y;
        if (refset_index == 1)
            y_orig_crlb = y_orig;
        end
        y  = y / max(y);  % normalize y-axis
        lw = 1.5;         % linewidth
        if refset_index == length(all_refs)
            lw = 2; 
        end
        plot(x,y,'Color',colors{refset_index},'LineWidth',lw);
        hold on; 
        grid on;

        % plot histogram for CRLB distribution
        if refset_index == 1

            % generate normalized histogram data
            [counts, edges] = histcounts(Cm,nbins);
            counts = counts / max(y_orig);

            % styling parameters
            bar(edges(1:end-1),counts,'FaceColor','#FEAE00','FaceAlpha',0.8,'EdgeAlpha',0);
            axis square;

            % axes labels
            title('\rm Glu','FontSize',24,'FontName','Calibri');
            xlabel('C_m (mM)');
            ylabel('Frequency');
            ax = gca;
            ax.FontSize = 18;
            ax.FontName = 'Calibri';
            ax.FontName = 'Calibri';
            set(gca,'yticklabel',[])
            set(gca,'ytick',[])
            ylim([0 1]);

            % set axis limit based on file
            xlim([10.9 11.8]);  % x-axis limits,
                                % dataset: 01_highFWHM_M_HC  [10.9 11.8]
                                % dataset: 02_lowFWHM_F_RRMS [14.0 15]
        end

        % calculate CV_adjusted
        if refset_index == length(all_refs) || refset_index == 1
            CV_normal   = std(Cm)/mean(Cm)
            CV_adjusted = sqrt(exp(std(log(Cm))^2)-1)
        end

        % refresh figure, pause for user input
        drawnow;
    end
    
    % save the figure as PNG
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMS.filename,'_glu_param.png'));
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMS.filename,'_glu_param.fig'),'fig');

    % --- show full width distribution --- 
    fig = figure;                      % initialize new figure
    fig.Position = [100,100,500,400];  % position, dimensions
    set(fig,'Color',[1 1 1]);          % set white figure background

    % distribution with CRLB of Am only
    [counts, edges] = histcounts(Cm_crlb,nbins);
    counts = counts / max(y_orig_crlb);
    bar( edges(1:end-1),counts, ...
         'FaceColor','#FEAE00', ...
         'FaceAlpha',0.9,       ...
         'EdgeAlpha',0);
    hold on; 

    % distribution of full propagated uncertainty
    [counts, edges] = histcounts(Cm,nbins);
    counts = counts / max(y_orig);
    bar( edges(1:end-1),counts, ...
         'FaceColor','#FEAE00', ...
         'FaceAlpha',0.4,       ...
         'EdgeAlpha',0);

    % style figure
    title('\rm Glu','FontSize',24,'FontName','Calibri'); 
    axis square;
    xlabel('C_m (mM)');
    ylabel('Frequency');
    ax = gca;
    ax.FontSize = 18; 
    ax.FontName = 'Calibri';
    ax.FontName = 'Calibri';

    % set x-axis limits based on dataset
    xlim([9.8 16]);
    ylim([0 1]);

    % determine best fit lognormal line
    h2 = histfit(Cm,nbins,'lognormal');
    h2(1).Visible = 'off';
    h2(2).Visible = 'off';
    x_lognormal = h2(2).XData; % density curve
    y_lognormal = h2(2).YData;
    y_lognormal = y_lognormal / max(y_lognormal); % normalize y-axis

    % plot best fit line
    plot(x,y,'Color','#FEAE00','LineWidth',2);
    plot(x_lognormal, y_lognormal, '--', 'Color', 'k', 'Linewidth', 1); % (lognormal)

    % save the figure as PNG
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMS.filename,'_glu_full.png'));
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMS.filename,'_glu_full.fig'),'fig');
end
