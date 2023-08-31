% ------------------------------------------------------------
% Script:  CmErrorCaseStudyMM.m
% Summary: Runs a Monte Carlo analysis to estimate the propagated
%          uncertainty of absolute quantification of healthy
%          volunteers from the referenced study below. 
%
% Inputs:  Quantification parameters are organized by subject, and
%          can be found here: 
%          ./data/mm_config/*
%
% Usage:   >> CmErrorCaseStudyMM.m
%
% Reference: 
%          Landheer K, Gajdošík M, Treacy M, Juchem C. Concentration 
%          and effective T2 relaxation times of macromolecules at 3T. 
%          Magn Reson Med. 2020;84(5):2327-2337. doi:10.1002/mrm.28282
% 
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
%
% Columbia University, 2023
% ------------------------------------------------------------
function CmErrorCaseStudyMM(config)
    close all; clc;
    disp('Running Case Study Analysis: MM dataset...');

    % analysis parameters
    n     = config.CmErrorCaseStudyMM.n;          % number of simulated Cm Values
    nbins = config.CmErrorCaseStudyMM.nbins;      % number of bins in simulated histogram
    ds    = CmUtils.load_dataset(config.CmErrorCaseStudyMM.filename);

    %% Perform Monte Carlo on macromolecule: MM2.04
    disp('simulating macromolecule: MM2.04...');

    % include all parameters in propagated uncertainty estimate
    % exclude f_csf since dependent on f_grey, f_white
    all_refs = {                             ...
        'Smm', 'Sw', 'T2mm','TE', 'TR',      ...
        'T1w_grey',  'T1w_white', 'T1w_csf', ...
        'T2w_grey',  'T2w_white', 'T2w_csf', ...
        'cw_grey',   'cw_white',  'cw_csf',  ...
        'f_grey' ,   'f_white' }; 
    curr_ref_set = {};

    % initialize figure
    fig = figure;
    fig.Position = [100,100,500,400]; % figure dimensions, size
    set(fig,'Color',[1 1 1]);         % set white figure background

    % simulate each subset of parameters
    for refset_index = 1:length(all_refs)

        % append to current parameter list
        curr_ref_set{end+1} = all_refs{refset_index};
        disp(strcat('Parameter of interest: ',curr_ref_set{end}));

        % load quantification parameters
        quant_params = {};
        quant_params.Smm       = ds.constants.macromolecules.mm2_04.Smm;
        quant_params.T2mm      = ds.constants.macromolecules.mm2_04.T2mm;
        quant_params.Sw        = ds.constants.water_ref.Sw;
        quant_params.T1w_grey  = ds.constants.water_ref.T1_grey;
        quant_params.T1w_white = ds.constants.water_ref.T1_white;
        quant_params.T1w_csf   = ds.constants.water_ref.T1_csf;
        quant_params.T2w_grey  = ds.constants.water_ref.T2_grey;
        quant_params.T2w_white = ds.constants.water_ref.T2_white;
        quant_params.T2w_csf   = ds.constants.water_ref.T2_csf;
        quant_params.cw_grey   = ds.constants.water_ref.cw_grey;
        quant_params.cw_white  = ds.constants.water_ref.cw_white;
        quant_params.cw_csf    = ds.constants.water_ref.cw_csf;
        quant_params.f_grey    = ds.constants.water_ref.f_grey;
        quant_params.f_white   = ds.constants.water_ref.f_white;
        quant_params.TR        = ds.constants.acquisition.TR;
        quant_params.TE        = ds.constants.acquisition.TE;

        % set quantification uncertaintites
        quant_uncert.dSmm       = ds.uncertainties.macromolecules.mm2_04.dSmm;
        quant_uncert.dT2mm      = ds.uncertainties.macromolecules.mm2_04.dT2mm;
        quant_uncert.dSw        = ds.uncertainties.water_ref.dSw;
        quant_uncert.dT1w_grey  = ds.uncertainties.water_ref.dT1_grey;
        quant_uncert.dT1w_white = ds.uncertainties.water_ref.dT1_white;
        quant_uncert.dT1w_csf   = ds.uncertainties.water_ref.dT1_csf;
        quant_uncert.dT2w_grey  = ds.uncertainties.water_ref.dT2_grey;
        quant_uncert.dT2w_white = ds.uncertainties.water_ref.dT2_white;
        quant_uncert.dT2w_csf   = ds.uncertainties.water_ref.dT2_csf;
        quant_uncert.dcw_grey   = ds.uncertainties.water_ref.dcw_grey;
        quant_uncert.dcw_white  = ds.uncertainties.water_ref.dcw_white;
        quant_uncert.dcw_csf    = ds.uncertainties.water_ref.dcw_csf;
        quant_uncert.df_grey    = ds.uncertainties.water_ref.df_grey;
        quant_uncert.df_white   = ds.uncertainties.water_ref.df_white;
        quant_uncert.df_csf     = ds.uncertainties.water_ref.df_csf;
        quant_uncert.dTR        = ds.uncertainties.acquisition.dTR;
        quant_uncert.dTE        = ds.uncertainties.acquisition.dTE;

        % set Monte Carlo simulation parameters
        sim_params = quant_params;      % handle constants with no variance in simulation

        % simulate a distribution for each reference parameter
        for ref_index = 1:length(curr_ref_set)

            % set current parameter
            curr_ref = curr_ref_set{ref_index};

            % set lower bound of truncation based on param
            if ismember(curr_ref,{'T2mm','T1w_grey','T1w_white','T1w_csf','T2w_grey','T2w_white','T2w_csf','TE','TR'})
                truncate_lower = 0.001; % min time of 1 ms
            else
                truncate_lower = 0; % default to non-negative A.U. values
            end

            % set upper bound for f_i, cannot have fraction > 1
            if ismember(curr_ref,{'f_grey'})
                if any(ismember(curr_ref_set,{'f_white'}))
                    truncate_upper = 1;
                else
                    truncate_upper = 1-quant_params.f_white;
                end
            else
                truncate_upper = inf;
            end

            % special case: f_white
            if ismember(curr_ref,{'f_white'})
                mu    = quant_params.(curr_ref);
                sigma = quant_uncert.(strcat('d',curr_ref));
                pd    = makedist('normal','mu',mu,'sigma',sigma);
                pd_random = zeros(1,n);
                for rand_index = 1:n
                    if mod(rand_index,20000) == 0
                        disp(strcat('rand_index:',num2str(rand_index)));
                    end
                    pd_truncated = truncate(pd,0,1-sim_params.('f_grey')(rand_index));
                    pd_random(rand_index) = random(pd_truncated);
                end
            else
                % simulate parameter distribution
                mu    = quant_params.(curr_ref);
                sigma = quant_uncert.(strcat('d',curr_ref));
                pd    = makedist('normal','mu',mu,'sigma',sigma);
                pd_truncated = truncate(pd,truncate_lower,truncate_upper); % truncate distribution
                pd_random    = random(pd_truncated,1,n); % generate random values
            end

            % replace simulation parameter constant with simulated distribution
            sim_params.(curr_ref) = pd_random;

        end

        % set value of f_csf based on f_grey and f_white
        quant_params.f_csf = 1 - (sim_params.('f_grey') + sim_params.('f_white'));
        quant_params.f_csf(quant_params.f_csf < 0) = 0; % avoid negative values
        sim_params.f_csf = quant_params.f_csf;

        % simulate absolute quantification (in M, since Cw is in M)
        Cm = CmUtils.absolute_quantification_mm(sim_params);
        Cm = Cm .* 1000; % convert to mM
        if (refset_index == 1)
            Cm_crlb = Cm;
        end

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
        y = y / max(y); % normalize y-axis
        if refset_index == length(all_refs)
            current_color = '#FEAE00';
        else
            current_color = '#5E5E5E';
        end
        switch config.CmErrorCaseStudyMM.filename
            case 'MM_01_B.json'
                if ismember(curr_ref,{'Smm','T2mm','f_grey','f_white'})
                    plot( x,y,                   ...
                          'Color',current_color, ...
                          'LineWidth',1.5); 
                    hold on; 
                end
            case 'MM_02_A.json'
                if ismember(curr_ref,{'Smm','T2mm','f_grey','f_white'})
                    plot( x,y,                   ...
                          'Color',current_color, ...
                          'LineWidth',1.5); 
                    hold on; 
                end
        end
        grid on;

        % plot histogram only for CRLB distribution
        if refset_index == 1

            % generate normalized histogram data
            [counts, edges] = histcounts(Cm,nbins);
            counts = counts / max(y_orig);

            % styling parameters
            bar( edges(1:end-1),counts, ...
                 'FaceColor','#FEAE00', ...
                 'FaceAlpha',0.8,       ...
                 'EdgeAlpha',0 );
            axis square;

            % axes labels
            title('\rm M2.04',   ...
                  'FontSize',24, ...
                  'FontName','Calibri');
            xlabel('C_m (mM)');
            ylabel('Frequency');
            ax = gca;
            ax.FontSize = 18;
            ax.FontName = 'Calibri';
            set(gca,'yticklabel',[])
            set(gca,'ytick',[])

            % set axis limits based on dataset
            switch config.CmErrorCaseStudyMM.filename
                case 'MM_01_B.json'
                    xlim([70.5 100]);
                    ylim([0 1]);
                case 'MM_02_A.json'
                    xlim([91.5 120]);
                    ylim([0 1]);
            end
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
    exportgraphics(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_M2.04_param.png'),'Resolution',2000);
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_M2.04_param.fig'),'fig');

    % --- show full width distribution --- 
    % initialize figure
    fig = figure;
    fig.Position = [100,100,500,400]; % position, dimensions
    set(fig,'Color',[1 1 1]);         % set white figure background

    % distribution with CRLB of Sm only
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
    bar(edges(1:end-1),counts, ...
        'FaceColor','#FEAE00', ...
        'FaceAlpha',0.4,       ...
        'EdgeAlpha',0);

    % style figure
    title('\rm M2.04', ...
        'FontSize',24, ...
        'FontName','Calibri'); 
    axis square;
    xlabel('C_m (mM)');
    ylabel('Frequency');
    ax = gca;
    ax.FontSize = 18;
    ax.FontName = 'Calibri';
    ylim([0 1]);
    xlim([40 140]);

    % determine best fit lognormal line
    h2 = histfit(Cm,nbins,'lognormal');
    h2(1).Visible = 'off';
    h2(2).Visible = 'off';
    x_lognormal = h2(2).XData;                    % density curve
    y_lognormal = h2(2).YData;
    y_lognormal = y_lognormal / max(y_lognormal); % normalize y-axis

    % plot best fit line 
    plot(x,y,'Color','#FEAE00','LineWidth',2); hold on;                 % (kernel)
    plot(x_lognormal, y_lognormal, '--', 'Color', 'k', 'Linewidth', 1); % (lognormal)

    % save the figure as PNG
    exportgraphics(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_M2.04_full.png'),'Resolution',2000);
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_M2.04_full.fig'),'fig');

    %% Perform Monte Carlo on macromolecule: MM0.92
    disp('simulating macromolecule: MM0.92...');

    % include all parameters in propagated uncertainty estimate
    % exclude f_csf since dependent on f_grey, f_white
    all_refs = {                            ...
        'Smm', 'Sw', 'T2mm', 'TE', 'TR',    ...
        'T1w_grey', 'T1w_white', 'T1w_csf', ...
        'T2w_grey', 'T2w_white', 'T2w_csf', ...
        'cw_grey',  'cw_white',  'cw_csf',  ...
        'f_grey' ,  'f_white'}; 
    curr_ref_set = {};

    % initialize figure
    fig = figure;
    fig.Position = [100,100,500,400]; % position, dimensions 
    set(fig,'Color',[1 1 1]);         % set white figure background

    for refset_index = 1:length(all_refs)

        % append to current parameter list
        curr_ref_set{end+1} = all_refs{refset_index};
        disp(strcat('Parameter of interest: ',curr_ref_set{end}));

        % load quantification parameters
        quant_params = {};
        quant_params.Smm       = ds.constants.macromolecules.mm0_92.Smm;
        quant_params.T2mm      = ds.constants.macromolecules.mm0_92.T2mm;
        quant_params.Sw        = ds.constants.water_ref.Sw;
        quant_params.T1w_grey  = ds.constants.water_ref.T1_grey;
        quant_params.T1w_white = ds.constants.water_ref.T1_white;
        quant_params.T1w_csf   = ds.constants.water_ref.T1_csf;
        quant_params.T2w_grey  = ds.constants.water_ref.T2_grey;
        quant_params.T2w_white = ds.constants.water_ref.T2_white;
        quant_params.T2w_csf   = ds.constants.water_ref.T2_csf;
        quant_params.cw_grey   = ds.constants.water_ref.cw_grey;
        quant_params.cw_white  = ds.constants.water_ref.cw_white;
        quant_params.cw_csf    = ds.constants.water_ref.cw_csf;
        quant_params.f_grey    = ds.constants.water_ref.f_grey;
        quant_params.f_white   = ds.constants.water_ref.f_white;
        quant_params.TR        = ds.constants.acquisition.TR;
        quant_params.TE        = ds.constants.acquisition.TE;

        % set quantification uncertaintites
        quant_uncert.dSmm       = ds.uncertainties.macromolecules.mm0_92.dSmm;
        quant_uncert.dT2mm      = ds.uncertainties.macromolecules.mm0_92.dT2mm;
        quant_uncert.dSw        = ds.uncertainties.water_ref.dSw;
        quant_uncert.dT1w_grey  = ds.uncertainties.water_ref.dT1_grey;
        quant_uncert.dT1w_white = ds.uncertainties.water_ref.dT1_white;
        quant_uncert.dT1w_csf   = ds.uncertainties.water_ref.dT1_csf;
        quant_uncert.dT2w_grey  = ds.uncertainties.water_ref.dT2_grey;
        quant_uncert.dT2w_white = ds.uncertainties.water_ref.dT2_white;
        quant_uncert.dT2w_csf   = ds.uncertainties.water_ref.dT2_csf;
        quant_uncert.dcw_grey   = ds.uncertainties.water_ref.dcw_grey;
        quant_uncert.dcw_white  = ds.uncertainties.water_ref.dcw_white;
        quant_uncert.dcw_csf    = ds.uncertainties.water_ref.dcw_csf;
        quant_uncert.df_grey    = ds.uncertainties.water_ref.df_grey;
        quant_uncert.df_white   = ds.uncertainties.water_ref.df_white;
        quant_uncert.df_csf     = ds.uncertainties.water_ref.df_csf;
        quant_uncert.dTR        = ds.uncertainties.acquisition.dTR;
        quant_uncert.dTE        = ds.uncertainties.acquisition.dTE;

        % set Monte Carlo simulation parameters
        sim_params = quant_params;      % handle constants with no variance in simulation

        % simulate a distribution for each reference parameter
        for ref_index = 1:length(curr_ref_set)

            % set current parameter
            curr_ref = curr_ref_set{ref_index};

            % set lower bound of truncation based on param
            if ismember(curr_ref,{'T2mm','T1w_grey','T1w_white','T1w_csf','T2w_grey','T2w_white','T2w_csf','TE','TR'})
                truncate_lower = 0.001; % min time of 1 ms
            else
                truncate_lower = 0; % default to non-negative A.U. values
            end

            % set upper bound for f_i, cannot have fraction > 1
            if ismember(curr_ref,{'f_grey'})
                if any(ismember(curr_ref_set,{'f_white'}))
                    truncate_upper = 1;
                else
                    truncate_upper = 1-quant_params.f_white;
                end
            else
                truncate_upper = inf;
            end

            % special case: f_white
            if ismember(curr_ref,{'f_white'})
                mu    = quant_params.(curr_ref);
                sigma = quant_uncert.(strcat('d',curr_ref));
                pd    = makedist('normal','mu',mu,'sigma',sigma);
                pd_random = zeros(1,n);
                for rand_index = 1:n
                    if mod(rand_index,20000) == 0
                        disp(strcat('rand_index:',num2str(rand_index)));
                    end
                    pd_truncated = truncate(pd,0,1-sim_params.('f_grey')(rand_index));
                    pd_random(rand_index) = random(pd_truncated);
                end
            else
                % simulate parameter distribution
                mu    = quant_params.(curr_ref);
                sigma = quant_uncert.(strcat('d',curr_ref));
                pd    = makedist('normal','mu',mu,'sigma',sigma);
                pd_truncated = truncate(pd,truncate_lower,truncate_upper); % truncate distribution
                pd_random = random(pd_truncated,1,n); % generate random values
            end

            % replace simulation parameter constant with simulated distribution
            sim_params.(curr_ref) = pd_random;

        end

        % set value of f_csf based on f_grey and f_white
        quant_params.f_csf = 1 - (sim_params.('f_grey') + sim_params.('f_white'));
        quant_params.f_csf(quant_params.f_csf < 0) = 0; % avoid negative values
        sim_params.f_csf = quant_params.f_csf;

        % simulate absolute quantification (in M, since Cw is in M)
        Cm = CmUtils.absolute_quantification_mm(sim_params);
        Cm = Cm .* 1000; % convert to mM
        if (refset_index == 1)
            Cm_crlb = Cm;
        end

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
        y = y / max(y); % normalize y-axis
        if refset_index == length(all_refs)
            current_color = '#0076BA';
        else
            current_color = '#5E5E5E';
        end
        if ismember(curr_ref,{'Smm','T2mm','cw_grey','f_grey','f_white'})
            plot(x,y,'Color',current_color,'LineWidth',1.5); hold on; 
        end
        grid on;

        % plot histogram only for CRLB distribution
        if refset_index == 1

            % generate normalized histogram data
            [counts, edges] = histcounts(Cm,nbins);
            counts = counts / max(y_orig);

            % styling parameters
            bar( edges(1:end-1),counts, ...
                 'FaceColor','#0076BA', ...
                 'FaceAlpha',0.8,       ...
                 'EdgeAlpha',0);
            axis square;

            % axes labels
            title('\rm M0.92','FontSize',24,'FontName','Calibri');
            xlabel('C_m (mM)');
            ylabel('Frequency');
            ax = gca;
            ax.FontSize = 18;
            ax.FontName = 'Calibri';
            set(gca,'yticklabel',[])
            set(gca,'ytick',[])

            % set axis limits based on dataset
            switch config.CmErrorCaseStudyMM.filename
                case 'MM_01_B.json'
                    xlim([25.09 31]);
                    ylim([0 1]);
                case 'MM_02_A.json'
                    xlim([28.15 35]);
                    ylim([0 1]);
            end
        end

        % calculate CV_adjusted
        if refset_index == length(all_refs) || refset_index == 1
            CV_normal   = std(Cm)/mean(Cm)
            CV_adjusted = sqrt(exp(std(log(Cm))^2)-1)
        end
        drawnow; % refresh figure, pause for user input
    end

    % save the figure as PNG
    exportgraphics(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_M0.92_param.png'),'Resolution',2000);
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_M0.92_param.fig'),'fig');

    % --- show full width distribution --- 
    % initialize figure
    fig = figure; fig.Position = [100,100,500,400];
    set(fig,'Color',[1 1 1]); % set white figure background

    % distribution with CRLB of Sm only
    [counts, edges] = histcounts(Cm_crlb,nbins);
    counts = counts / max(y_orig_crlb);
    bar(edges(1:end-1),counts,'FaceColor','#0076BA','FaceAlpha',0.9,'EdgeAlpha',0);
    hold on; 

    % distribution of full propagated uncertainty
    [counts, edges] = histcounts(Cm,nbins);
    counts = counts / max(y_orig);
    bar(edges(1:end-1),counts,'FaceColor','#0076BA','FaceAlpha',0.4,'EdgeAlpha',0);

    % style figure
    title('\rm M0.92',   ...
          'FontSize',24, ...
          'FontName','Calibri'); 
    axis square;
    xlabel('C_m (mM)');
    ylabel('Frequency');
    ax = gca;
    ax.FontSize = 18;
    ax.FontName = 'Calibri';
    ylim([0 1]);
    xlim([20 38]);

    % determine best fit lognormal line
    h2 = histfit(Cm,nbins,'lognormal');
    h2(1).Visible = 'off';
    h2(2).Visible = 'off';
    x_lognormal   = h2(2).XData; % density curve
    y_lognormal   = h2(2).YData;
    y_lognormal   = y_lognormal / max(y_lognormal); % normalize y-axis

    % plot best fit line
    plot(x,y,'Color','#0076BA','LineWidth',2);                          % (kernel)
    plot(x_lognormal, y_lognormal, '--', 'Color', 'k', 'Linewidth', 1); % (lognormal)

    % save the figure as PNG
    exportgraphics(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_M0.92_full.png'),'Resolution',2000);
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_M0.92_full.fig'),'fig');


    %% Perform Monte Carlo on metabolite: tNAA
    disp('simulating macromolecule: tNAA...');

    % include all parameters in propagated uncertainty estimate
    % exclude f_csf since dependent on f_grey, f_white
    all_refs = { 
        'Sm', 'Sw', 'T1m', 'T2m', 'TE', 'TR', ...
        'T1w_grey', 'T1w_white', 'T1w_csf',   ...
        'T2w_grey', 'T2w_white', 'T2w_csf',   ...
        'cw_grey',  'cw_white',  'cw_csf',    ...
        'f_grey' ,  'f_white'};
    curr_ref_set = {};

    % initialize figure
    fig = figure;
    fig.Position = [100,100,500,400]; % position, dimensions
    set(fig,'Color',[1 1 1]);         % set white figure background

    for refset_index = 1:length(all_refs)

        % append to current parameter list
        curr_ref_set{end+1} = all_refs{refset_index};
        disp(strcat('Parameter of interest: ',curr_ref_set{end}));

        % load quantification parameters
        quant_params = {};
        quant_params.Sm        = ds.constants.metabolites.tNAA.Sm;
        quant_params.T1m       = ds.constants.metabolites.tNAA.T1m;
        quant_params.T2m       = ds.constants.metabolites.tNAA.T2m;
        quant_params.Sw        = ds.constants.water_ref.Sw;
        quant_params.T1w_grey  = ds.constants.water_ref.T1_grey;
        quant_params.T1w_white = ds.constants.water_ref.T1_white;
        quant_params.T1w_csf   = ds.constants.water_ref.T1_csf;
        quant_params.T2w_grey  = ds.constants.water_ref.T2_grey;
        quant_params.T2w_white = ds.constants.water_ref.T2_white;
        quant_params.T2w_csf   = ds.constants.water_ref.T2_csf;
        quant_params.cw_grey   = ds.constants.water_ref.cw_grey;
        quant_params.cw_white  = ds.constants.water_ref.cw_white;
        quant_params.cw_csf    = ds.constants.water_ref.cw_csf;
        quant_params.f_grey    = ds.constants.water_ref.f_grey;
        quant_params.f_white   = ds.constants.water_ref.f_white;
        quant_params.TR        = ds.constants.acquisition.TR;
        quant_params.TE        = ds.constants.acquisition.TE;

        % set quantification uncertaintites
        quant_uncert.dSm        = ds.uncertainties.metabolites.tNAA.dSm;
        quant_uncert.dT1m       = ds.uncertainties.metabolites.tNAA.dT1m;
        quant_uncert.dT2m       = ds.uncertainties.metabolites.tNAA.dT2m;
        quant_uncert.dSw        = ds.uncertainties.water_ref.dSw;
        quant_uncert.dT1w_grey  = ds.uncertainties.water_ref.dT1_grey;
        quant_uncert.dT1w_white = ds.uncertainties.water_ref.dT1_white;
        quant_uncert.dT1w_csf   = ds.uncertainties.water_ref.dT1_csf;
        quant_uncert.dT2w_grey  = ds.uncertainties.water_ref.dT2_grey;
        quant_uncert.dT2w_white = ds.uncertainties.water_ref.dT2_white;
        quant_uncert.dT2w_csf   = ds.uncertainties.water_ref.dT2_csf;
        quant_uncert.dcw_grey   = ds.uncertainties.water_ref.dcw_grey;
        quant_uncert.dcw_white  = ds.uncertainties.water_ref.dcw_white;
        quant_uncert.dcw_csf    = ds.uncertainties.water_ref.dcw_csf;
        quant_uncert.df_grey    = ds.uncertainties.water_ref.df_grey;
        quant_uncert.df_white   = ds.uncertainties.water_ref.df_white;
        quant_uncert.df_csf     = ds.uncertainties.water_ref.df_csf;
        quant_uncert.dTR        = ds.uncertainties.acquisition.dTR;
        quant_uncert.dTE        = ds.uncertainties.acquisition.dTE;

        % set Monte Carlo simulation parameters
        sim_params = quant_params;      % handle constants with no variance in simulation

        % simulate a distribution for each reference parameter
        for ref_index = 1:length(curr_ref_set)

            % set current parameter
            curr_ref = curr_ref_set{ref_index};

            % set lower bound of truncation based on param
            if ismember(curr_ref,{'T1m', 'T2m','T1w_grey','T1w_white','T1w_csf','T2w_grey','T2w_white','T2w_csf','TE','TR'})
                truncate_lower = 0.001; % min time of 1 ms
            else
                truncate_lower = 0; % default to non-negative A.U. values
            end

            % set upper bound for f_i, cannot have fraction > 1
            if ismember(curr_ref,{'f_grey'})
                if any(ismember(curr_ref_set,{'f_white'}))
                    truncate_upper = 1;
                else
                    truncate_upper = 1-quant_params.f_white;
                end
            else
                truncate_upper = inf;
            end

            % special case: f_white
            if ismember(curr_ref,{'f_white'})
                mu    = quant_params.(curr_ref);
                sigma = quant_uncert.(strcat('d',curr_ref));
                pd    = makedist('normal','mu',mu,'sigma',sigma);
                pd_random = zeros(1,n);
                for rand_index = 1:n
                    if mod(rand_index,20000) == 0
                        disp(strcat('rand_index:',num2str(rand_index)));
                    end
                    pd_truncated = truncate(pd,0,1-sim_params.('f_grey')(rand_index));
                    pd_random(rand_index) = random(pd_truncated);
                end
            else
                % simulate parameter distribution
                mu    = quant_params.(curr_ref);
                sigma = quant_uncert.(strcat('d',curr_ref));
                pd    = makedist('normal','mu',mu,'sigma',sigma);
                pd_truncated = truncate(pd,truncate_lower,truncate_upper); % truncate distribution
                pd_random    = random(pd_truncated,1,n); % generate random values
            end

            % replace simulation parameter constant with simulated distribution
            sim_params.(curr_ref) = pd_random;
        end

        % set value of f_csf based on f_grey and f_white
        quant_params.f_csf = 1 - (sim_params.('f_grey') + sim_params.('f_white'));
        quant_params.f_csf(quant_params.f_csf < 0) = 0; % avoid negative values
        sim_params.f_csf = quant_params.f_csf;

        % simulate absolute quantification (in M, since Cw is in M)
        Cm = CmUtils.absolute_quantification_water(sim_params);
        Cm = Cm .* 1000; % convert to mM
        if (refset_index == 1)
            Cm_crlb = Cm;
        end

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
        y = y / max(y); % normalize y-axis
        if refset_index == length(all_refs)
            current_color = '#0076BA';
        else
            current_color = '#5E5E5E';
        end
        if ismember(curr_ref,{'Sm','T1m','T2m','cw_grey','f_grey','f_white'})
            plot(x,y,'Color',current_color,'LineWidth',1.5);
            hold on; 
        end
        grid on;

        % plot histogram only for CRLB distribution
        if refset_index == 1

            % generate normalized histogram data
            [counts, edges] = histcounts(Cm,nbins);
            counts = counts / max(y_orig);

            % styling parameters
            bar( edges(1:end-1),counts, ...
                 'FaceColor','#0076BA', ...
                 'FaceAlpha',0.8,       ...
                 'EdgeAlpha',0);
            axis square;

            % axes labels
            title('\rm tNAA','FontSize',24,'FontName','Calibri');
            xlabel('C_m (mM)');
            ylabel('Frequency');
            ax = gca;
            ax.FontSize = 18;
            ax.FontName = 'Calibri'; 
            set(gca,'yticklabel',[])
            set(gca,'ytick',[])

            % set axis limits based on dataset
            switch config.CmErrorCaseStudyMM.filename
                case 'MM_01_B.json'
                    xlim([11.125 14]);
                    ylim([0 1]);
                case 'MM_02_A.json'
                    xlim([10.68 13.5]);
                    ylim([0 1]);
            end
        end

        % calculate CV_adjusted
        if refset_index == length(all_refs) || refset_index == 1
            CV_normal   = std(Cm)/mean(Cm)
            CV_adjusted = sqrt(exp(std(log(Cm))^2)-1)
        end
        drawnow; % refresh figure
    end

    % save the figure as PNG
    exportgraphics(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_tNAA_param.png'),'Resolution',2000);
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_tNAA_param.fig'),'fig');

    % --- show full width distribution --- 
    fig = figure;                      % initialize figure
    fig.Position = [100,100,500,400];  % position, dimensions
    set(fig,'Color',[1 1 1]);          % set white figure background

    % distribution with CRLB of Sm only
    [counts, edges] = histcounts(Cm_crlb,nbins);
    counts = counts / max(y_orig_crlb);
    bar(edges(1:end-1),counts,'FaceColor','#0076BA','FaceAlpha',0.9,'EdgeAlpha',0);
    hold on; 

    % distribution of full propagated uncertainty
    [counts, edges] = histcounts(Cm,nbins);
    counts = counts / max(y_orig);
    bar(edges(1:end-1),counts,'FaceColor','#0076BA','FaceAlpha',0.4,'EdgeAlpha',0);

    % style figure
    title( '\rm tNAA',    ...
           'FontSize',24, ...
           'FontName','Calibri'); 
    axis square;
    xlabel('C_m (mM)');
    ylabel('Frequency');
    ax = gca;
    ax.FontSize = 18;
    ax.FontName = 'Calibri';
    ylim([0 1]);
    xlim([8 16]);

    % determine best fit lognormal line
    h2 = histfit(Cm,nbins,'lognormal');
    h2(1).Visible = 'off';
    h2(2).Visible = 'off';
    x_lognormal   = h2(2).XData; % density curve
    y_lognormal   = h2(2).YData;
    y_lognormal   = y_lognormal / max(y_lognormal); % normalize y-axis

    % plot best fit line
    plot(x,y,'Color','#0076BA','LineWidth',2);
    plot(x_lognormal, y_lognormal, '--', 'Color', 'k', 'Linewidth', 1); % (lognormal)

    % save the figure as PNG
    exportgraphics(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_tNAA_full.png'),'Resolution',2000);
    saveas(fig,strcat(config.paths.res_dir,'mc_case_study_', ...
        config.CmErrorCaseStudyMM.filename,'_tNAA_full.fig'),'fig');
end
