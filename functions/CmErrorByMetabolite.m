% ------------------------------------------------------------
% Function: CmErrorByMetabolite
% Summary:  returns a simulated concentration error, dCm from MR spectroscopy
%           absolute quantification, Cm for a specific metabolite. This
%           includes the full propagated error estimate using reported
%           standard deviations of all absolute quantification parameters,
%           as listed below. Both the concentration estimate Cm and
%           propagated error, dCm are returned and organized by metabolite.
%
%           This implementation allows users to specify any of the following 
%           absolute quantification parameters, including:
%               Am    - the LCM quantification amplitude
%               Sw    - signal water amplitude
%               T1m   - metabolite longitudinal relaxation time (seconds)
%               T2m   - metabolite transverse relaxation time (seconds)
%               TE    - sequence echo time (seconds)
%               TR    - sequence repetition time (seconds)
%               T1_i  - longitudinal relaxation time of water in grey matter,
%                       white matter and CSF
%               T2_i  - transerver relaxation time of water in grey matter,
%                       white matter and CSF
%               f_i   - fractional estimate of grey and white matter
%                       in the interrogated voxel of interest
%               c_wi  - fractional concentration of water in grey matter,
%                       white matter, and CSF
%
% Inputs:   params - struct, includes shared and metabolite-specific
%               constants and uncertainties, according to the following
%               data model: 
%           {
%               "options": {
%                   "metabolites":   {"naa","glu","gln", ...},
%                   "includedTerms": {"Am", "Sw", "T1m", ...}
%               },
%               "constants": {
%                   "params": {"Sw":1, "TE":0.02, ...},
%                   "metabolites": {
%                       "naa": {"Am": 0.0001, "T1m": 1.38, "T2m": 0.295},
%                       "glu": {"Am": 0.0001, "T1m": 1.27, "T2m": 0.122}
%                       ...
%                   },
%               },
%               "errors": {
%                   "upper_bound": {
%                       "params": {"dSw":0.05, ...},
%                       "metabolites": {
%                           "naa": { "dAm": ##, "dT1m": ##, "dT2m": ## },
%                           "glu": { "dAm": ##, "dT1m": ##, "dT2m": ## },
%                       }
%                   },
%                   "lower_bound": {...}
%               },
%               "covariances": {
%                   "upper_bound": { "Am-Sw": ##, ... },
%                   "lower_bound": { "Am-Sw": ##, ... }
%               }
%           }
% 
% Outputs:  Cm  - struct, estimated absolute concentration by metabolite
%           dCm - struct, estimated propagated error by metabolite
%           md  - struct, metadata containing results & parameters from
%                 intermediary steps.
%
% Notes:    This function uses the equation the absolute quantification found
%           in the following publication:
%
%           [1] K. Landheer, M. Gajdošík, and C. Juchem, “A semi‐LASER, 
%           single‐voxel spectroscopic sequence with a minimal echo time of 
%           20.1 ms in the human brain at 3 T,” NMR in Biomedicine, vol. 33, 
%           no. 9, Sep. 2020, doi: 10.1002/nbm.4324.
%
%           Covariances are in the form x-y, which is equivalent to y-x,
%           since covar matrix is symmetric, e.g. Am-Sw == Sw-Am. System
%           defaults a covariance value to 0 is unspecified.
%
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
%
% Columbia University, 2023
% ------------------------------------------------------------
function [Cm, dCm, md] = CmErrorByMetabolite(params)
    
    % exit on invalid params
    required = {'options.metabolites', 'options.includedTerms', ...
        'errors.params', 'errors.metabolites', ...
        'constants.params', 'constants.metabolites'};
    for i = 1:length(required)
        splitf = split(required{i},'.');
        checkfield = params;
        for j = 1:length(splitf)
            currfield = splitf{j};
            if ~isfield(checkfield,splitf{j})
                strcat({'Field: '},splitf{j},' is not contained in params, exiting...')
                Cm = 0; dCm = 0; return; 
            end
            checkfield = checkfield.(currfield);
        end
    end
    
    % save parameters as metadata
    md.params = params;

    % toggle metabolites to include/exclude from error estimate
    METABOLITES = {'cho','cr','naa','gaba','gln','glu'};
    for i = 1:length(METABOLITES)
        if isfield(params.options,'metabolites') && ...
            any(strcmp(params.options.metabolites,METABOLITES{i}))
            meta_opts.(METABOLITES{i}) = true;  % include metabolite
        else
            meta_opts.(METABOLITES{i}) = false; % exclude metabolite
        end
    end
    
    % toggle included terms (e.g. Am, T1m, T2m)
    TERMS = { 'Am','T1m','T2m','Sw','TE','TR', ...
              'f_csf','f_grey','f_white',      ...
              'cw_csf','cw_grey','cw_white',   ...
              'T1_csf', 'T1_grey', 'T1_white', ...
              'T2_csf', 'T2_grey', 'T2_white'
    };
    for i = 1:length(TERMS)
        if (isfield(params.options,'includedTerms') && ...
            any(strcmp(params.options.includedTerms,TERMS{i})))
            term_opts.(TERMS{i}) = true;  % include term in dCm
        else
            term_opts.(TERMS{i}) = false; % exclude term from dCm
        end
    end

    % set constants, errors
    c = params.constants.params;
    er = params.errors.params;
    c_meta = params.constants.metabolites;
    er_meta = params.errors.metabolites;
    
    % load covariances
    covariances = params.covariances;
    
    % calculate dCm error for each metabolite
    metas = fieldnames(meta_opts);
    for i = 1:length(metas)
        
        % ensure metabolite fields exist
        if meta_opts.(metas{i}) == true && isfield(c_meta,metas{i}) ...
            && isfield(er_meta,metas{i})
            disp(strcat('calculating quantification error for metabolite:',metas{i}));

            % use CmPartials class to calculate partials
            curr_consts = c;
            curr_consts.Am  = c_meta.(metas{i}).Am;
            curr_consts.T1m = c_meta.(metas{i}).T1m;
            curr_consts.T2m = c_meta.(metas{i}).T2m;
            cmp_obj = CmPartials(curr_consts);
            partials = cmp_obj.calc_partials.partials;
            Cm.(metas{i}) = cmp_obj.Cm;
            
            % save partials as metadata
            md.results.partials.(metas{i}) = partials;

            % calculate Gaussian error propagation, dCM
            % using Jacobian and Covariance matrix calculation 
            dCm_meta = []; 
            for ref_i = 1:numel(er_meta.(metas{i}).dAm)

                % display progress in console
                if mod(ref_i,100) == 0
                    disp(strcat('reference Am index:',int2str(ref_i),'-of-',int2str(numel(er_meta.(metas{i}).dAm))));
                end

                % construct Jacobian matrix
                J = [];
                all_fields = fieldnames(term_opts);
                fields = []; % restrict to included fields
                for field_i = 1:numel(all_fields)
                    if term_opts.(all_fields{field_i}) == true
                        fields = [fields all_fields(field_i)];
                    end
                end
                for all_field_i = 1:numel(fields)
                    if strcmp(fields{all_field_i},'Am')
                        J = [J partials.(strcat(fields{all_field_i}))];
                    else
                        J = [J partials.(strcat(fields{all_field_i}))];
                    end
                end

                % construct covariance matrix
                SIGMA = [];
                for sigma_i = 1:numel(fields)
                    for sigma_j = 1:numel(fields)
                        if sigma_i == sigma_j
                            if strcmp(fields{sigma_i},'Am')
                                SIGMA(sigma_i,sigma_j) = er_meta.(metas{i}).(strcat('d',fields{sigma_i}))(ref_i)^2;
                            elseif strcmp(fields{sigma_i},'T1m') || strcmp(fields{sigma_i},'T2m')
                                SIGMA(sigma_i,sigma_j) = er_meta.(metas{i}).(strcat('d',fields{sigma_i}))^2;
                            else
                                SIGMA(sigma_i,sigma_j) = er.(strcat('d',fields{sigma_i}))^2;
                            end
                        else
                            % matrix is symmetric, so cov(x,y) is same as cov(y,x)
                            % use either value if defined in JSON config
                            cov_field_forward = strcat(fields{sigma_i},'_',fields{sigma_j});
                            cov_field_reverse = strcat(fields{sigma_j},'_',fields{sigma_i});
                            if isfield(covariances,cov_field_forward)
                                SIGMA(sigma_i,sigma_j) = covariances.(cov_field_forward);
                            elseif isfield(covariances,cov_field_reverse)
                                SIGMA(sigma_i,sigma_j) = covariances.(cov_field_reverse); 
                            else
                                % default to 0 covariance
                                SIGMA(sigma_i,sigma_j) = 0;
                            end
                        end
                    end
                end
                
                % save relevant metadata
                md.results.sigma.(metas{i}){ref_i} = SIGMA;
                md.results.J.(metas{i}){ref_i} = J;

                % calculate Gaussian error propagation, dCM,
                % based on Jacobian and Covariance matrix
                dCm_meta = [dCm_meta sqrt(J*SIGMA*J')];
            end
            dCm.(metas{i}) = dCm_meta;
        end
    end
    
    % save results as metadata
    md.results = Cm;
    md.results = dCm;
end