% ------------------------------------------------------------
% Script:  CmErrorMonteCarloByParam.m
% Summary: performs a Monte Carlo simulation to estimate the propagated
%          error of absolute quantification estimate (Cm) from MR
%          spectroscopy datasets in the human brain, isolating the effect
%          of each quantification parameter uncertainty on the overall
%          propagated error. A histogram of resulting Cm values is
%          generated for each parameter, displaying the distributions
%          resulting from the highest and lowest parameter uncertainties.
%
% Usage:   >> CmErrorMonteCarloByParam
%
% Inputs:  @file - ../data/CmErrorParams.json
%          contains parameter constants and reported standard deviations
%          from relevant MRS literature. fields such as "covariances" are
%          ignored in this script, since Monte Carlo simulations are not
%          dependent on partial derivatives of parameter covariances.
% 
% Output:  @file - a PNG image with the filename monte-carlo-[param].png
%          is rendered and saved to disk. 
%
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
%
% Columbia University, 2023
% ------------------------------------------------------------
function CmErrorMonteCarloByParam(config)
    close all; clc;
    disp('Running Monte Carlo error propagation analysis by parameter...');

    % script parameters
    n  = config.CmErrorMonteCarloByParam.n;              % number of simulated Cm values
    all_params = config.CmErrorMonteCarloByParam.params; % quantification parameters
    ds = CmUtils.load_dataset('CmErrorParams.json');     % load constants

    % colors for each parameter, shaded error ranges as HEX
    colors = {
        '#929292', '#0076BA', '#76D6FF', ...  % grey, blue, light blue,
        '#FEAE00', '#9830A0',            ...  % orange, purple,
        '#1DB100', '#1DB100', '#1DB100', ...  % green (x3),
        '#0076BA', '#0076BA', '#0076BA', ...  % blue (x3),
        '#76D6FF', '#76D6FF', '#76D6FF', ...  % light blue (x3),
        '#EE220C', '#EE220C', '#EE220C'  ...  % red (x3)
    }';

    % run analysis on each quantification parameter
    for i = 1:length(all_params)

        % set current parameter (e.g. Am, Sw)
        param = all_params{i};
        disp(strcat('Analyzing quantification parameter:',param));

        % set random variable, Normally distributed
        % ~N(parameter value, standard deviation)
        mu    = ds.constants.(param);
        sigma_upper = ds.errors.upper_bound.(strcat('d',param));
        sigma_lower = ds.errors.lower_bound.(strcat('d',param));
        r_upper = normrnd(mu,sigma_upper,[1 n]);
        r_lower = normrnd(mu,sigma_lower,[1 n]);

        % use truncated normalized random variable
        % for variables: T1m, T2m, T1_i, T2_i, f_i
        if strcmp(param,'T1m') || strcmp(param,'T1_white') || ...
            strcmp(param,'T1_grey')  || strcmp(param,'T1_csf')  || ...
            strcmp(param,'T2_white') || strcmp(param,'T2_grey') || ...
            strcmp(param,'T2_csf')
            pd_upper = makedist('normal','mu',mu,'sigma',sigma_upper);
            pd_lower = makedist('normal','mu',mu,'sigma',sigma_lower);
            t_upper  = truncate(pd_upper,1e-3,inf);
            t_lower  = truncate(pd_lower,1e-3,inf);
            r_upper  = random(t_upper,1,n);
            r_lower  = random(t_lower,1,n);
        end
        if strcmp(param,'T2m')
            pd_upper = makedist('normal','mu',mu,'sigma',sigma_upper);
            pd_lower = makedist('normal','mu',mu,'sigma',sigma_lower);
            t_upper  = truncate(pd_upper,1e-3,inf);
            t_lower  = truncate(pd_lower,1e-3,inf);
            r_upper  = random(t_upper,1,n);
            r_lower  = random(t_lower,1,n);
        end
        if strcmp(param,'f_csf') || strcmp(param,'f_grey') || strcmp(param,'f_white')
            pd_upper = makedist('normal','mu',mu,'sigma',sigma_upper);
            pd_lower = makedist('normal','mu',mu,'sigma',sigma_lower);
            t_upper  = truncate(pd_upper,0,1);
            t_lower  = truncate(pd_lower,0,1);
            r_upper  = random(t_upper,1,n);
            r_lower  = random(t_lower,1,n);
        end
        if strcmp(param,'cw_csf') || strcmp(param,'cw_grey') || strcmp(param,'cw_white')
            pd_upper = makedist('normal','mu',mu,'sigma',sigma_upper);
            pd_lower = makedist('normal','mu',mu,'sigma',sigma_lower);
            t_upper  = truncate(pd_upper,0,inf);
            t_lower  = truncate(pd_lower,0,inf);
            r_upper  = random(t_upper,1,n);
            r_lower  = random(t_lower,1,n);
        end

        % set parameter to randomly generated value
        params_upper = ds.constants;
        params_lower = ds.constants;
        params_upper.(param) = r_upper;
        params_lower.(param) = r_lower;

        % calculate Cm
        cmobj_upper = CmPartials(params_upper);
        cmobj_lower = CmPartials(params_lower);
        disp('calculating Cm using upper bounds...');
        Cm_upper = cmobj_upper.Cm .* 1000;    % convert to mM
        disp('calculating Cm using lower bounds...');
        Cm_lower = cmobj_lower.Cm .* 1000;
        disp('Done.');

        % check lognormal CV for comparison
        cv_raw_lognormal = sqrt(exp(std(log(Cm_upper))^2)-1)
        std_upper = std(Cm_upper)

        % plot histogram of variable of interest, calculated Cm values
        fig = figure;                      % initialize figure
        fig.Position = [100,100,500,400];  % size and position
        nbins = 1000;                      % histogram bins
        set(fig,'Color',[1 1 1]);          % set figure background to white
        [Cml_counts, edgesL] = histcounts(Cm_lower-mean(Cm_lower),nbins);
        [Cmu_counts, edgesU] = histcounts(Cm_upper-mean(Cm_upper),nbins);
        Cml_counts = Cml_counts / max(Cml_counts);
        Cmu_counts = Cmu_counts / max(Cmu_counts);
        bar(edgesL(1:end-1),Cml_counts, 'FaceColor',colors{i},...
            'FaceAlpha',0.8,'EdgeAlpha',0,'BarWidth', 1); hold on;
        bar(edgesU(1:end-1),Cmu_counts, 'FaceColor',colors{i},...
            'FaceAlpha',0.3,'EdgeAlpha',0,'BarWidth', 1);
        title(CmUtils.format_param_as_string(strcat('\rm',param)));
        axis square;
        xlabel('$C_m-\bar{C}_m$','Interpreter','Latex');
        ylabel('Frequency');
        ax = gca;
        ax.FontSize = 18;
        ax.FontName = 'Calibri';
        set(gca,'yticklabel',[])
        set(gca,'ytick',[])
        set(gca,'YLim',[0 1])

        % adjust axis limits per group
        xlim([-1.0 1.0]);
        if strcmp(param,'T2_csf') || strcmp(param,'T2_white')
            xlim([-.1 .1]);
        end
        if strcmp(param,'T2_grey') || strcmp(param,'T1_grey') || strcmp(param,'T1_white') || strcmp(param,'T1_csf')
            xlim([-.5 .5]);
        end
        if strcmp(param,'T1m')
            xlim([-0.68 .68]);
        end

        % add annotation of std/sample_mean as textbox
        dim = [.62 .78 .1 .1]; % location: upper right
        str = { strcat("$\widehat{CV}_{upper}=",num2str(compose("%.3f",std(double(Cm_upper))/mean(double(Cm_upper)))),"$"),...
                strcat("$\widehat{CV}_{lower}=",num2str(compose("%.3f",std(double(Cm_lower))/mean(double(Cm_lower)))),"$")};
        if (strcmp(param,'T2m'))
            str = {strcat("$\widehat{CV}_{lower}=",num2str(compose("%.3f",std(double(Cm_lower))/mean(double(Cm_lower)))),"$")};
        end
        a = annotation('textbox',dim,'interpreter','latex','String',str, ...
            'verticalalignment','middle','horizontalalignment','center');
        a.FontSize = 16;

        % save the figure
        saveas(fig,strcat(config.paths.res_dir,'monte-carlo-',param,'.png'));
    end
end