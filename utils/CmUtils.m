% ------------------------------------------------------------
% Class:   CmUtils
% Summary: A set of utility functions shared by analysis scripts performing
%          error propagation simulations and calculations. All
%          miscellaneous functions shared by scripts are implemented here
%          to prevent implementation redundancy. 
% 
% Usage:   myUtils = CmUtils()
%          myUtils.load_dataset('CmErrorMetabs.json');
%
% Inputs:  N/A
%
% Methods: (1) load_dataset - loads a JSON dataset and returns as a struct
%          (2) format_param_as_string - returns a formatted string for chart
%                   parameters, handling underscores and subscripts
%          (3) foramt_metab_as_string - returns a formatted string for
%                   chart metabolites, handling underscores and subscripts          
%
% Author:  Ronald Instrella, Columbia University
% PI:      Christoph Juchem, Ph.D.
% ------------------------------------------------------------
classdef CmUtils
    
    % public methods (endpoints)
    % no class contructor is necessary, define all as Static methods
    methods(Static) 
        
        % loads and returns a JSON dataset from local
        function [ds] = load_dataset(file_name)
            
            fid = fopen(file_name);
            raw = fread(fid,inf); 
            str = char(raw');
            fclose(fid); 
            ds = jsondecode(str);
            
        end
        
        % returns a formatted string for chart parameters, underscores & symbols
        function [param_str] = format_param_as_string(param)
            
            param_str = strrep(param,    'Am',       'A_{m}');
            param_str = strrep(param_str,'Sw',       'S_{w}');
            param_str = strrep(param_str,'T1m',      '{\it T}_{1m}');
            param_str = strrep(param_str,'T2m',      '{\it T}_{2m}');
            param_str = strrep(param_str,'f_csf',    'f_{csf}');
            param_str = strrep(param_str,'f_grey',   'f_{grey}');
            param_str = strrep(param_str,'f_white',  'f_{white}');
            param_str = strrep(param_str,'T1_csf',   '{\it T}_{1w,csf}');
            param_str = strrep(param_str,'T1_grey',  '{\it T}_{1w,grey}');
            param_str = strrep(param_str,'T1_white', '{\it T}_{1w,white}');
            param_str = strrep(param_str,'T2_csf',   '{\it T}_{2w,csf}');
            param_str = strrep(param_str,'T2_grey',  '{\it T}_{2w,grey}');
            param_str = strrep(param_str,'T2_white', '{\it T}_{2w,white}');
            param_str = strrep(param_str,'cw_csf',   'C_{w,csf}');
            param_str = strrep(param_str,'cw_grey',  'C_{w,grey}');
            param_str = strrep(param_str,'cw_white', 'C_{w,white}');
            param_str = strrep(param_str,'TE',       '{\it T}_{E}');
            param_str = strrep(param_str,'TR',       '{\it T}_{R}');
        end
        
        % returns a formatted string for chart metabolites
        function [param_str] = format_metab_as_string(param)
            
            param_str = strrep(param,    'gln', 'Gln');
            param_str = strrep(param_str,'glu', 'Glu');
            param_str = strrep(param_str,'cho', 'Cho');
            param_str = strrep(param_str,'gaba','GABA');
            param_str = strrep(param_str,'naa', 'NAA');
            param_str = strrep(param_str,'cr',  'Cr');
            
        end
        
        % returns the absolute concentration of generic metabolite with 
        % LCM amplitude Am; T1 and T2 corrections are applied,
        % quantification is performed using an internal tCr reference
        % input params: Am, T1m, T2m, S_tCr, T1_tCr, T2_tCr, TE, TR, C_tCr
        function [Cm] = absolute_quantification_tCr(params)

            % T1 correction
            T1_corr = (1-exp(-params.TR./params.T1_tCr)) ./ (1-exp(-params.TR./params.T1m));
            
            % T2 correction 
            T2_corr = exp(-params.TE./params.T2_tCr) ./ exp(-params.TE./params.T2m);
            
            % absolute concentration estimate of metabolite m
            Cm = (params.Am ./ params.S_tCr) .* T1_corr .* T2_corr .* params.C_tCr;
        end

        % returns the absolute concentration of generic macromolecule with
        % LCM amplitude Amm, using a tissue-specific water reference
        % T1 and T2 corrections are applied for water; 
        % T2 correction is applied for macromolecule
        % input params: 
        %    Amm, T2mm, Sw, 
        %    f_csf, f_grey, f_white,
        %    cw_grey, cw_white, cw_csf, 
        %    T1w_grey, T1w_white, T1w_csf, 
        %    T2w_grey, T2w_white, T2w_csf
        function [Cmm] = absolute_quantification_mm(params)
            
            % T1 water corrections
            T1w_grey_corr  = (1-exp(-params.TR ./ params.T1w_grey));
            T1w_white_corr = (1-exp(-params.TR ./ params.T1w_white));
            T1w_csf_corr   = (1-exp(-params.TR ./ params.T1w_csf));
            
            % T2 water corrections
            T2w_grey_corr  = exp(params.TE ./ params.T2w_grey);
            T2w_white_corr = exp(params.TE ./ params.T2w_white);
            T2w_csf_corr   = exp(params.TE ./ params.T2w_csf);
            
            % T2 macromolecule corrections
            T2mm_correction = exp(params.TE ./ params.T2mm);
            
            % Tissue specific expressions
            corr_grey  = params.f_grey  .* T1w_grey_corr  .* params.cw_grey  ./ T2w_grey_corr;
            corr_white = params.f_white .* T1w_white_corr .* params.cw_white ./ T2w_white_corr;
            corr_csf   = params.f_csf   .* T1w_csf_corr   .* params.cw_csf   ./ T2w_csf_corr;
            sum_corrections = corr_grey + corr_white + corr_csf;
            
            % absolute concentration estimate
            Cmm = (params.Amm ./ (params.Sw .* (1-params.f_csf))) .* T2mm_correction .* sum_corrections;
        end
        
        % returns the absolute concentration of generic metabolic m with
        % LCM amplitude Am, using tissue-specific water reference
        % T1 and T2 corrections are applied for water and the metabolite
        function [Cm] = absolute_quantification_water(params)
            
            % T1 water corrections
            T1w_grey_corr  = (1-exp(-params.TR ./ params.T1w_grey));
            T1w_white_corr = (1-exp(-params.TR ./ params.T1w_white));
            T1w_csf_corr   = (1-exp(-params.TR ./ params.T1w_csf));
            
            % T2 water corrections
            T2w_grey_corr  = exp(params.TE ./ params.T2w_grey);
            T2w_white_corr = exp(params.TE ./ params.T2w_white);
            T2w_csf_corr   = exp(params.TE ./ params.T2w_csf);
            
            % T1 and T2 metabolite corrections
            T1m_correction = 1 ./ (1-exp(-params.TR ./ params.T1m));
            T2m_correction = exp(params.TE ./ params.T2m);
            
            % Tissue specific expressions
            corr_grey  = params.f_grey  .* T1w_grey_corr  .* params.cw_grey  ./ T2w_grey_corr;
            corr_white = params.f_white .* T1w_white_corr .* params.cw_white ./ T2w_white_corr;
            corr_csf   = params.f_csf   .* T1w_csf_corr   .* params.cw_csf   ./ T2w_csf_corr;
            sum_corrections = corr_grey + corr_white + corr_csf;
            
            % absolute concentration estimate
            Cm = (params.Am ./ (params.Sw .* (1-params.f_csf))) .* T1m_correction .* T2m_correction .* sum_corrections;
        end
    end
end






