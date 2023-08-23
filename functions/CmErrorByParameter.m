% ------------------------------------------------------------
% Function: CmErrorByParameter
% Summary:  returns a simulated concentration error, dCm from MR spectroscopy
%           absolute quantification estimate, Cm. Options are provided to
%           either display the full propagated error, including covariances and
%           uncertainties of various system and sequence parameters, or view an 
%           approximation using a subset of all parameter uncertainties.
%
%           This implementation varies the uncertainty of a selected
%           parameter of interest (options.reference_var) over a specified
%           range of values (options.error_range), holding all other
%           parameters and uncertainties constant, and calculates the
%           estimated absolute quantification (Cm) and propagated error
%           (dCm) over the specified error range.
%
% Inputs:   params - struct, includes all constants, parameter
%               uncertainties, and methods options, according to the
%               following data model: 
%           {
%               "options": {
%                   "reference_var": 'Am',
%                   "error_range":   (0:.01:1)
%               },
%               "constants":   {Am: 0.0001, Sw: 1, ...},
%               "errors":      {dAm: 0.05, dSw: 0.06, ...},
%               "covariances": {f_csf-f_grey: 0.001, ...}
%               "description": "this is a test dataset"
%           }
%
% Outputs:  Cm  - vector, estimated absolute concentration
%           dCm - vector, estimated propagated error
%           md  - struct, metadata containing data from intermediary steps
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
%           The following sequence & system parameters are supported: 
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
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
%
% Columbia University, 2023
% ------------------------------------------------------------
function [Cm, dCm, md] = CmErrorByParameter(params)

    % exit on invalid params
    required = {'options', 'errors', 'constants'};
    for i = 1:length(required)
        if ~isfield(params,required{i})
            Cm = 0; dCm = 0; return; 
        end
    end
    
    % save relevant metadata
    md.params = params;

    % error frame of refernce (x-axis)
    options.reference_var = params.options.reference_var;
    options.error_range = params.options.error_range;

    % define terms to include/exclude from error estimate
    TERMS = { ...
        'Am','Sw','T1m','T2m','TE','TR', ...
        'f_grey', 'f_white', 'f_csf',    ...
        'cw_csf', 'cw_grey', 'cw_white', ...
        'T1_csf', 'T1_grey', 'T1_white', ...
        'T2_csf', 'T2_grey', 'T2_white'
    };
    for i = 1:length(TERMS)
        if ~isfield(params.options,'includedTerms')
            term_opts.(TERMS{i}) = true;  % default behavior
        elseif any(strcmp(params.options.includedTerms,TERMS{i}))
            term_opts.(TERMS{i}) = true;  % include term
        else
            term_opts.(TERMS{i}) = false; % exclude term
        end
    end

    % get constants, absolute errors
    c  = params.constants; 
    er = params.errors;
    covariances = params.covariances;

    % replace primary error with range
    reference_err      = strcat('d',options.reference_var);
    er.(reference_err) = options.error_range;

    % create instance of CmPartials class
    cmp_obj = CmPartials(c);

    % calculate total concentration estimate, Cm
    Cm = (c.Am*exp(c.TE/c.T2m)*cmp_obj.Dm)/(c.Sw*(1-c.f_csf))*cmp_obj.f_sum;
    
    % calculate partial derivatives
    partials = cmp_obj.calc_partials.partials;
    
    % save partials as metadata
    md.results.Cm = Cm;
    md.results.partials = partials;
    
    % calculate dCm for each reference value
    dCm = [];
    for ref_i = 1:numel(er.(reference_err))
    
        % construct Jacobian matrix
        J = [];
        all_fields = fieldnames(term_opts);
        fields = []; % restrict to included fields
        for field_i = 1:numel(all_fields)
            if term_opts.(all_fields{field_i}) == true
                fields = [fields all_fields(field_i)];
            end
        end
        for i = 1:numel(fields)
            J = [J partials.(fields{i})];
        end

        % construct covariance matrix
        SIGMA = [];
        for i = 1:numel(fields)
            for j = 1:numel(fields)
                if i == j
                    if strcmp(fields{i},options.reference_var)
                        index_to_assign = ref_i;
                    else
                        index_to_assign = 1;
                    end
                    SIGMA(i,j) = er.(strcat('d',fields{i}))(index_to_assign)^2;
                else
                    % matrix is symmetric, so cov(x,y) is same as cov(y,x)
                    % use either value if defined in JSON config
                    cov_field_forward = strcat(fields{i},'-',fields{j});
                    cov_field_reverse = strcat(fields{j},'-',fields{i});
                    if isfield(covariances,cov_field_forward)
                        SIGMA(i,j) = covariances.(cov_field_forward);
                    elseif isfield(covariances,cov_field_reverse)
                        SIGMA(i,j) = covariances.(cov_field_reverse); 
                    else
                        % default to 0 covariance
                        SIGMA(i,j) = 0;
                    end
                end
            end
        end
        
        % save relevant metadata
        md.results.sigma(:,:,ref_i) = SIGMA;
        
        % calculate Gaussian error propagation, dCM,
        % based on Jacobian and Covariance matrix
        dCm = [dCm sqrt(J*SIGMA*J')];
    end
    
    % save relevant metadata
    md.results.dCm = dCm;
end