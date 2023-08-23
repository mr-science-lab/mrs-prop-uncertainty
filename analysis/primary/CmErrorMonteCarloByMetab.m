% ------------------------------------------------------------
% Script:  CmErrorMonteCarloByMetab.m
% Summary: performs a Monte Carlo simultation to estimate the propagated
%          error of absolute quantification estimate (Cm) from MR spectroscopy
%          datasets in the human brain, for a user-specified metabolite.
%          A histogram of resulting Cm values is generated for each metabolite, 
%          displaying the distributions resulting from the highest and lowest
%          set of reported parameter uncertainties.
%
% Usage:   >> CmErrorMonteCarloByMetab
%
% Inputs:  @file - ../data/CmErrorMetabs.json
%          contains parameter constants and reported standard deviations
%          from relevant MRS literature. fields such as "covariances"
%          are ignored in this script, since the Monte Carlo simulations are
%          not dependent on partial derviatives or parameter covariances.
%
% Output:  @file - a PNG image with filename monte-carlo-[metabolite].png
%          is rendered and saved to disk.
%
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
%
% Columbia University, 2023
% ------------------------------------------------------------
function CmErrorMonteCarloByMetab(config)
    close all; clc;
    disp('Running Monte Carlo error propagation analysis by metabolite...');

    % script parameters
    n      = config.CmErrorMonteCarloByMetab.n;            % number of simulated Cm values
    nbins  = config.CmErrorMonteCarloByMetab.nbins;        % number of histogram bins
    ds     = CmUtils.load_dataset('CmErrorMetabs.json');   % load constants
    metas  = config.CmErrorMonteCarloByMetab.metabolites;  % metabolites to simulate
    colors = config.CmErrorMonteCarloByMetab.colors;       % colors for each metabolite plot

    % include all parameters in propagated error estimate
    % note: order of f_i is relevant, exclude f_csf
    refs = config.CmErrorMonteCarloByMetab.parameters;

    % run analysis on each metabolite
    for metab_index = 1:length(metas)

        % get current metabolite, initialize constants
        disp(strcat('simulating metabolite:',metas{metab_index}));
        curr_metab = metas{metab_index};
        params_upper = {}; params_lower = {};

        % generate a random variable for each parameter;
        % need to iterate instead of vector operations (as in Monte Carlo
        % analysis by parameter), due to the dependency between random
        % variables f_csf, f_white, and f_grey. 
        for ref_index = 1:length(refs)

            % set current parameter
            curr_ref = refs{ref_index};

            % set random variable for each metabolite-specific parameter
            if ismember(curr_ref,{'Am','T1m','T2m'})

                % case: metabolite-specific parameters
                mu = ds.metabolites.(curr_metab).(curr_ref);
                sigma_upper = ds.errors.upper_bound.metabolites.(curr_metab).(strcat('d',curr_ref));
                sigma_lower = ds.errors.lower_bound.metabolites.(curr_metab).(strcat('d',curr_ref));

                % hold Am constant between upper and lower
                % use the upper bound for Am across simulations
                if ismember(curr_ref,{'Am'})
                    sigma_lower = ds.errors.upper_bound.metabolites.(curr_metab).('dAm');
                end

                % truncate T1m, T2m at 0, avoid negative numbers
                if ismember(curr_ref,{'T1m','T2m'})
                    pd_upper = makedist('normal','mu',mu,'sigma',sigma_upper);
                    pd_lower = makedist('normal','mu',mu,'sigma',sigma_lower);
                    t_upper  = truncate(pd_upper,1e-2,inf);
                    t_lower  = truncate(pd_lower,1e-2,inf);
                    r_upper  = random(t_upper,1,n);
                    r_lower  = random(t_lower,1,n);
                else % case: Am
                    r_upper = normrnd(mu,sigma_upper,[1 n]);
                    r_lower = normrnd(mu,sigma_lower,[1 n]);
                end

            % set random variable for parameters shared between metabolites 
            elseif ismember(curr_ref,{'Sw','TE','TR','f_csf','f_grey','f_white', ...
                'cw_csf', 'cw_grey', 'cw_white','T1_csf', 'T1_grey', 'T1_white', ...
                'T2_csf', 'T2_grey', 'T2_white'})

                % case: shared constants
                mu = ds.constants.(curr_ref);
                sigma_upper = ds.errors.upper_bound.params.(strcat('d',curr_ref));
                sigma_lower = ds.errors.lower_bound.params.(strcat('d',curr_ref));

                % restrict f_i values accordingly
                if ismember(curr_ref,{'f_grey'})

                    % set f_grey to truncated normal random variable
                    % restricted to domain [0,1]
                    pd_upper = makedist('normal','mu',mu,'sigma',sigma_upper);
                    pd_lower = makedist('normal','mu',mu,'sigma',sigma_lower);
                    t_upper  = truncate(pd_upper,0,1);
                    t_lower  = truncate(pd_lower,0,1);
                    r_upper  = random(t_upper,1,n);
                    r_lower  = random(t_lower,1,n);

                elseif ismember(curr_ref,{'f_white'})

                    % set f_white to truncate normal random variable
                    % restricted to domain [0,1-f_grey]
                    disp('generating random values of f_white...');
                    pd_upper = makedist('normal','mu',mu,'sigma',sigma_upper);
                    pd_lower = makedist('normal','mu',mu,'sigma',sigma_lower);
                    r_upper = zeros(1,n);
                    r_lower = zeros(1,n);
                    for rand_index = 1:n
                        if mod(rand_index,100000) == 0
                            disp(strcat('rand_index:',num2str(rand_index)));
                        end
                        t_upper = truncate(pd_upper,0,1-params_upper.('f_grey')(rand_index));
                        t_lower = truncate(pd_lower,0,1-params_lower.('f_grey')(rand_index));
                        r_upper(rand_index) = random(t_upper);
                        r_lower(rand_index) = random(t_lower);
                    end
                    disp('done.'); % takes time, so print to console when complete.

                else

                    % truncate Sw, TE, TR, T1wi, T2wi, cwi to a truncated
                    % normal random variable from [0 inf], must be non-negative. 
                    pd_upper = makedist('normal','mu',mu,'sigma',sigma_upper);
                    pd_lower = makedist('normal','mu',mu,'sigma',sigma_lower);
                    t_upper  = truncate(pd_upper,0,inf);
                    t_lower  = truncate(pd_lower,0,inf);
                    r_upper  = random(t_upper,1,n);
                    r_lower  = random(t_lower,1,n);

                end
            end

            % set param to random variable
            params_upper.(curr_ref) = r_upper;
            params_lower.(curr_ref) = r_lower;
        end

        % Generate distribution using only CRLB
        curr_ref = 'Am';
        mu = ds.metabolites.(curr_metab).(curr_ref);
        sigma_crlb = ds.errors.upper_bound.metabolites.(curr_metab).(strcat('d',curr_ref));
        r_crlb = normrnd(mu,sigma_crlb,[1 n]);
        params_crlb = ds.constants;
        params_crlb.T1m = ds.metabolites.(curr_metab).T1m;
        params_crlb.T2m = ds.metabolites.(curr_metab).T2m;
        params_crlb.(curr_ref) = r_crlb;
        cmobj_crlb = CmPartials(params_crlb);
        Cm_crlb = cmobj_crlb.Cm;
        Cm_crlb = Cm_crlb .* 1000;

        % calcluate Cm using simulated params
        cmobj_upper = CmPartials(params_upper);
        cmobj_lower = CmPartials(params_lower);
        Cm_upper = cmobj_upper.Cm;      % in Molar (c_wi is in Molar)
        Cm_lower = cmobj_lower.Cm; 
        Cm_upper = Cm_upper .* 1000;    % convert to mM
        Cm_lower = Cm_lower .* 1000;

        % plot histogram of calculated Cm values
        fig = figure;
        fig.Position = [100,100,500,400]; 
        set(fig,'Color',[1 1 1]); % set white figure background

        % normalize by max amplitude of each distribution 
        [Cml_counts, edgesL] = histcounts(Cm_lower-mean(Cm_crlb),nbins);
        [Cmu_counts, edgesU] = histcounts(Cm_upper-mean(Cm_crlb),nbins);
        [Cmc_counts, edgesC] = histcounts(Cm_crlb-mean(Cm_crlb), nbins);
        Cml_counts = Cml_counts / max(Cml_counts);
        Cmu_counts = Cmu_counts / max(Cmu_counts);
        Cmc_counts = Cmc_counts / max(Cmc_counts);

        % render histograms
        bar(edgesL(1:end-1),Cml_counts, 'FaceColor',colors{metab_index},...
            'FaceAlpha',0.8,'EdgeAlpha',0,'BarWidth', 1); hold on;
        bar(edgesU(1:end-1),Cmu_counts, 'FaceColor',colors{metab_index},...
            'FaceAlpha',0.4,'EdgeAlpha',0,'BarWidth', 1); 
        bar(edgesC(1:end-1),Cmc_counts, 'FaceColor','black',...
            'FaceAlpha',0.4,'EdgeAlpha',0,'BarWidth', 1); 

        % format title, axes and fonts
        title(CmUtils.format_metab_as_string(strcat('\rm',curr_metab)), ...
            'FontSize',24,'FontName','Calibri'); 
        axis square;
        xlabel('$C_m-\bar{C}_m$','Interpreter','Latex');
        ylabel('Frequency');
        ax = gca;
        ax.FontSize = 18;
        ax.FontName = 'Calibri';
        ax.FontName = 'Calibri';
        set(gca,'YLim',[0 1], 'YTick',[0:1:1])
        set(gca,'XLim',[-5 5],'XTick',[-5:2.5:5])
        ytickformat('%1.f');

        % use for annotation of CV on each plot
        dim = [.62 .78 .1 .1];  % set position, size of textbox
        disp(strcat('Metabolite: ',metas{metab_index}));
        CV_upper = sqrt(exp(std(log(Cm_upper))^2)-1)
        CV_lower = sqrt(exp(std(log(Cm_lower))^2)-1)
        CV_crlb  = sqrt(exp(std(log(Cm_crlb))^2)-1)
        disp(strcat('Mean of CRLB distribution: ',num2str(mean(Cm_crlb))))
        disp(strcat('Mean of lower distribution:',num2str(mean(Cm_lower))))
        disp(strcat('Mean of upper distribution:',num2str(mean(Cm_upper))))

        % include standard CV for gaussian distributions
        CV_upper_normal = std(Cm_upper)/mean(Cm_upper)
        CV_lower_normal = std(Cm_lower)/mean(Cm_lower)

        % use lognormal definition of coefficient of variation
        str = { strcat("$\widehat{CV}_{high}=",num2str(compose("%.3f",std(Cm_upper)/mean(Cm_upper))),"$"),...
                strcat("$\widehat{CV}_{low}=",num2str(compose("%.3f",std(Cm_lower)/mean(Cm_lower))),"$")};
        a.FontSize  = 16;

        % save the figure as PNG
        saveas(fig,strcat(config.paths.res_dir,'monte-carlo-',curr_metab,'.png'));
    end
end