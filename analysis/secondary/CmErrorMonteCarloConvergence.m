% ------------------------------------------------------------
% Script:  CmErrorMonteCarloConvergence.m
% Summary: Plots the coefficient of variation as a function of number of
%          Monte Carlo iterations, for 6 metabolites. 
%          Displayed value should stabilize at higher iterations.
%
% Inputs:  None. Simulation parameters are found here:
%          ./data/CmErrorParams.json
%
% Usage:   >> CmErrorMonteCarloConvergence.m
% 
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
%
% Columbia University, 2023
% ------------------------------------------------------------
function CmErrorMonteCarloConvergence(config)
    close all; clc;

    % set script parameters
    meta     = config.CmErrorMonteCarloConvergence.metabolite; % selected metabolite
    meta_str = config.CmErrorMonteCarloConvergence.meta_str;   % metabolite string labels
    n        = [1:2:1000];                                     % vector for total number of iterations

    % track results
    cvu = [];  % coefficient of variation, upper bound
    cvl = [];  % coefficient of variation, lower bound
    cvc = [];  % coefficient of variation, CRLB

    % run for each number of iterations
    for i = [1:length(n)]
        ni = n(i);
        [cvu_curr, cvl_curr, cvc_curr] = run_cm_monte_carlo_by_metabolite(meta,ni);
        cvu(i) = cvu_curr;      % upper bound
        cvl(i) = cvl_curr;      % lower bound
        cvc(i) = cvc_curr;      % CRLB
    end

    % plot convergence curve
    fig = figure; 
    fig.Position = [100,100,500,400];                           % figure dimensions & position 
    set(fig,'Color',[1 1 1]);                                   % set background color to white
    plot(n,cvu,'k-','LineWidth',2);                             % upper bound (black)
    hold on; grid on;                                           % include all 3 plots on single axes
    plot(n,cvl,'-','Color','#0076BA','LineWidth',2);  hold on;  % lower bound (blue)
    plot(n,cvc,'-','Color',[.8 .8 .8],'LineWidth',2); hold on;  % CRLB (grey)
    legend('upper bound', 'lower bound', 'CRLB (ref)');         % add figure legend
    title(meta_str,'FontWeight','normal','Fontsize',24);        % label by metabolite
    ax          = gca; 
    ax.FontSize = 20; 
    ax.FontName = 'Calibri';
    xl          = xlabel(strcat('Iteration'));
    yl          = ylabel('CV, \deltaC_m/C_m [%]');
    set(xl,'FontSize',20);
    set(yl,'FontSize',20);

    % save to disk
    saveas(fig,strcat(config.paths.res_dir,'convergence-by-metabolite-',meta,'.fig'));
    saveas(fig,strcat(config.paths.res_dir,'convergence-by-metabolite-',meta,'.png'));

    
    % ------------------------------------------------------------
    % FUNCTION: run_cm_monte_carlo_by_metabolite
    % SUMMARY:  generates the coefficent of variation for the upper bound, 
    %           lower bound and CRLB monte carlo simulations at CRLB = 0, 
    %           for the number of simulations, n
    % ------------------------------------------------------------
    function [CV_upper, CV_lower, CV_crlb] = run_cm_monte_carlo_by_metabolite(meta,n)

        disp(['running simulation n = ' num2str(n)]);
        ds = CmUtils.load_dataset('CmErrorMetabs.json');  % load constants
        metas = {meta};

        % include all parameters in propagated error estimate
        % order of f_i is relevant, exclude f_csf
        refs  = { ...
            'Am','T1m','T2m','Sw','TE','TR', ...
            'f_grey', 'f_white',             ... 
            'cw_csf', 'cw_grey', 'cw_white', ...
            'T1_csf', 'T1_grey', 'T1_white', ...
            'T2_csf', 'T2_grey', 'T2_white'
        };

        % run analysis on each metabolite 
        for metab_index = 1:length(metas)

            % get current metabolite, initialize constants
            disp(strcat('simulating metabolite:',metas{metab_index}));
            curr_metab = metas{metab_index};
            params_upper = {}; params_lower = {};

            % generate a random variable for each parameter;
            % need to iterate instead of vector operations (as in Monte Carlo
            % analysis by parameter), due to the dependency between random
            % variables f_white, and f_grey. 
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
                        t_upper  = truncate(pd_upper,9e-3,inf);
                        t_lower  = truncate(pd_lower,9e-3,inf);
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
                            if mod(rand_index,20000) == 0
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

            % add annotation
            disp(strcat('Metabolite: ',metas{metab_index}));
            CV_upper = sqrt(exp(std(log(Cm_upper))^2)-1)
            CV_lower = sqrt(exp(std(log(Cm_lower))^2)-1)
            CV_crlb  = sqrt(exp(std(log(Cm_crlb))^2)-1)
            disp(strcat('Mean of CRLB distribution: ',num2str(mean(Cm_crlb))))
            disp(strcat('Mean of lower distribution:',num2str(mean(Cm_lower))))
            disp(strcat('Mean of upper distribution:',num2str(mean(Cm_upper))))
        end
    end
end