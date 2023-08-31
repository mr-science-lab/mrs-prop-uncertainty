% ------------------------------------------------------------
% Class:   CmPartials
% Summary: Provides a series of utilities and static methods to calculate
%          partial derivatives of terms used in the equation for absolute
%          quantification in MR Spectroscopy, which converts metabolite 
%          concentrations from arbitrary to absolute units (e.g. mmol).
%          All derviatives dCm/dxi are calculated using analytical solutions.
% 
% Usage:   myObj = CmPartials(constants)
%          myObj.calc_partials
% 
% Notes:   @constants - struct, must contain the following fields with
%           float or int: Sm, Sw, TE, TR ,T1m, T2m, f_grey, f_white, f_csf,
%           T1_grey, T1_white, T1_csf, T2_grey, T2_white, T2_csf,
%           cw_grey, cw_white, cw_grey.
%          
%          This Class uses the equation for absolute quantification found
%          in the following publication:
%
%          [1] K. Landheer, M. Gajdošík, and C. Juchem, “A semi‐LASER, 
%          single‐voxel spectroscopic sequence with a minimal echo time of 
%          20.1 ms in the human brain at 3 T,” NMR in Biomedicine, vol. 33, 
%          no. 9, Sep. 2020, doi: 10.1002/nbm.4324.
%
% Author:  Ronald Instrella, Columbia University
% PI:      Christoph Juchem, Ph.D.
% ------------------------------------------------------------
classdef CmPartials
    
    % class properties
    properties
        % Includes constants for calculating partials 
        consts = struct; 
    end
    
    % public methods 
    % all methods assume scalar values in current implementation
    methods
        
        % class constructor
        % @myConsts - sets quantification constants
        function obj = CmPartials(myConsts)
            % set instance variable
            obj.consts = myConsts;
        end
        
        % returns struct of partials
        function res = calc_partials(obj)
            
            % calls all public partial functions below
            partials.Sm       = obj.partial_Sm;
            partials.Sw       = obj.partial_Sw;
            partials.T1m      = obj.partial_T1m;
            partials.T2m      = obj.partial_T2m;
            partials.TE       = obj.partial_TE;
            partials.TR       = obj.partial_TR;
            partials.f_grey   = obj.partial_fgrey;
            partials.f_white  = obj.partial_fwhite;
            partials.f_csf    = obj.partial_fcsf;
            partials.T1_grey  = obj.partial_T1grey;
            partials.T1_white = obj.partial_T1white;
            partials.T1_csf   = obj.partial_T1csf;
            partials.T2_grey  = obj.partial_T2grey;
            partials.T2_white = obj.partial_T2white;
            partials.T2_csf   = obj.partial_T2csf;
            partials.cw_grey  = obj.partial_cwgrey;
            partials.cw_white = obj.partial_cwwhite;
            partials.cw_csf   = obj.partial_cwcsf;
            
            % return under struct field 'partials'
            res.partials = partials;
        end
        
        % dCm/dSm
        function val = partial_Sm(obj)
            
            val = (exp(obj.consts.TE/obj.consts.T2m) * obj.Dm) / ...
                  (obj.consts.Sw*(obj.consts.f_white + obj.consts.f_grey)) * obj.f_sum;
        end
        
        % dCm/dSw
        function val = partial_Sw(obj)
            
            val = (-1/(obj.consts.Sw^2)) * ...
                  (obj.consts.Sm * exp(obj.consts.TE / obj.consts.T2m) * ...
                  obj.Dm) / (obj.consts.f_white + obj.consts.f_grey) * obj.f_sum;
        end
        
        % dCm/dT1m
        function val = partial_T1m(obj)
            
            % group terms as scalar constant
            R_dT1m = (obj.consts.Sm * exp(obj.consts.TE/obj.consts.T2m)) ...
                     / (obj.consts.Sw * (obj.consts.f_white + obj.consts.f_grey)) * obj.f_sum;
            val    = R_dT1m * (1-exp(-obj.consts.TR/obj.consts.T1m))^-2 ...
                     * (exp(-obj.consts.TR/obj.consts.T1m) ...
                     * (obj.consts.TR/(obj.consts.T1m^2)));
        end
        
        % dCm/dT2m
        function val = partial_T2m(obj)
            
            val = exp(obj.consts.TE / obj.consts.T2m)   * ...
                  (-obj.consts.TE / (obj.consts.T2m^2)) * ...
                  (obj.consts.Sm*obj.Dm)/(obj.consts.Sw * ...
                  (obj.consts.f_white + obj.consts.f_grey)) * obj.f_sum;
        end
        
        % dCm/df_grey
        function val = partial_fgrey(obj)
            R   = obj.shared_fi_terms;
            val = R.fshared * ( obj.consts.f_white * R.fgrey ...
                  - obj.consts.f_white * R.fwhite - R.fcsf);
        end
        
        %dCm/df_white
        function val = partial_fwhite(obj)
            R   = obj.shared_fi_terms;
            val = R.fshared * (obj.consts.f_grey * R.fwhite ...
                  - obj.consts.f_grey * R.fgrey  - R.fcsf);
        end
        
        %dCm/df_csf
        function val = partial_fcsf(obj)
            R   = obj.shared_fi_terms;
            val = R.fshared * (R.fcsf - obj.consts.f_white ...
                * R.fgrey - obj.consts.f_grey * R.fwhite);
        end
        
        % dCm/dTE
        function val = partial_TE(obj)
            
            % f_i terms
            f_csf_dTE   = (obj.consts.f_csf   * obj.consts.cw_csf   / obj.Dw_csf   );
            f_grey_dTE  = (obj.consts.f_grey  * obj.consts.cw_grey  / obj.Dw_grey  );
            f_white_dTE = (obj.consts.f_white * obj.consts.cw_white / obj.Dw_white );
            
            % group scalar factor
            R_dTE = (obj.consts.Sm*obj.Dm)/(obj.consts.Sw*(obj.consts.f_white + obj.consts.f_grey));
            val   = R_dTE * ( exp( obj.consts.TE/obj.consts.T2m ...
                    - obj.consts.TE / obj.consts.T2_csf ) ...
                    * (1/obj.consts.T2m - 1/obj.consts.T2_csf)  * f_csf_dTE  ...
                    + exp(obj.consts.TE/obj.consts.T2m-obj.consts.TE/obj.consts.T2_grey )  ...
                    * (1/obj.consts.T2m - 1/obj.consts.T2_grey) * f_grey_dTE ...
                    + exp(obj.consts.TE/obj.consts.T2m-obj.consts.TE/obj.consts.T2_white ) ...
                    * (1/obj.consts.T2m - 1/obj.consts.T2_white)* f_white_dTE );
        end
        
        % dCm/dTR
        function val = partial_TR(obj)
            
            % shared terms
            R_dTR = obj.consts.Sm*exp(obj.consts.TE/obj.consts.T2m) ...
                    /(obj.consts.Sw*(1-obj.consts.f_csf));
            Mz = (1-exp(-obj.consts.TR/obj.consts.T1m));
            
            % dTR grey matter terms
            f_grey_dTR = (obj.Dm*exp(-obj.consts.TR/obj.consts.T1_grey) ...
                       / (obj.consts.T1_grey)) ...
                       + ((-1)*(exp(-obj.consts.TR/obj.consts.T1m)...
                       / ( obj.Dw_grey * Mz^2 * obj.consts.T1m)));
            f_grey_dTR_scaled = f_grey_dTR * obj.consts.f_grey * obj.consts.cw_grey ...
                              / exp(obj.consts.TE/obj.consts.T2_grey);
            
            % dTR white matter terms
            f_white_dTR = (obj.Dm*exp(-obj.consts.TR/obj.consts.T1_white) ...
                        /(obj.consts.T1_white)) ...
                        + ((-1)*(exp(-obj.consts.TR/obj.consts.T1m) ...
                        /(obj.Dw_white * Mz^2 * obj.consts.T1m)));
            f_white_dTR_scaled = f_white_dTR * obj.consts.f_white * obj.consts.cw_white ...
                               / exp(obj.consts.TE/obj.consts.T2_white);
            
            % dTR CSF terms
            f_csf_dTR = (obj.Dm*exp(-obj.consts.TR/obj.consts.T1_csf) ...
                        /(obj.consts.T1_csf)) ...
                        + ((-1)*(exp(-obj.consts.TR/obj.consts.T1m) ...
                        /(obj.Dw_csf * Mz^2*obj.consts.T1m)));
            f_csf_dTR_scaled = f_csf_dTR * obj.consts.f_csf * obj.consts.cw_csf ...
                             / exp(obj.consts.TE/obj.consts.T2_csf);
            
            % final expression of dCm/dTR
            val = R_dTR * (f_csf_dTR_scaled + f_grey_dTR_scaled + f_white_dTR_scaled);
        end
        
        % dCm/dT1_grey
        function val = partial_T1grey(obj)

            R_dT1_grey = (obj.consts.Sm*exp(obj.consts.TE/obj.consts.T2m) ...
                       * obj.Dm * obj.consts.f_grey * obj.consts.cw_grey) ...
                       / (obj.consts.Sw * (obj.consts.f_white + obj.consts.f_grey) ...
                       * exp(obj.consts.TE / obj.consts.T2_grey));
                   
            val = -((obj.consts.TR * exp(-obj.consts.TR / obj.consts.T1_grey)) ...
                / (obj.consts.T1_grey^2) ) * R_dT1_grey;
        end
        
        % dCm/dT1_white
        function val = partial_T1white(obj)

            R_dT1_white = (obj.consts.Sm*exp(obj.consts.TE/obj.consts.T2m) ...
                        * obj.Dm*obj.consts.f_white*obj.consts.cw_white) ...
                        / (obj.consts.Sw*(obj.consts.f_white + obj.consts.f_grey) ...
                        * exp(obj.consts.TE/obj.consts.T2_white));
            
            val = -((obj.consts.TR * exp(-obj.consts.TR/obj.consts.T1_white)) ... 
                / (obj.consts.T1_white^2) ) * R_dT1_white;
        end
        
        % dCm/dT1_csf
        function val = partial_T1csf(obj)
            
            R_dT1_csf = (obj.consts.Sm*exp(obj.consts.TE/obj.consts.T2m) ...
                      * obj.Dm*obj.consts.f_csf*obj.consts.cw_csf) ...
                      / (obj.consts.Sw*(obj.consts.f_white + obj.consts.f_grey) ...
                      * exp(obj.consts.TE/obj.consts.T2_csf));
          
            val = -((obj.consts.TR * exp(-obj.consts.TR/obj.consts.T1_csf)) ...
                / (obj.consts.T1_csf^2) ) * R_dT1_csf;
        end
        
        % dCm/dT2_grey
        function val = partial_T2grey(obj)
            
            % all T2w_i exprssions should have (1-f_csf) in denominator
            R_dT2_grey  = (obj.consts.Sm*exp(obj.consts.TE/obj.consts.T2m) ...
                        * obj.Dm * obj.consts.f_grey ...
                        * obj.consts.cw_grey) ...
                        / (obj.consts.Sw * (obj.consts.f_white + obj.consts.f_grey) * obj.Dw_grey);
                    
            val = ((obj.consts.TE * exp(-obj.consts.TE/obj.consts.T2_grey)) ...
                / (obj.consts.T2_grey^2)) * R_dT2_grey;
        end
        
        % dCm/dT2_white
        function val = partial_T2white(obj)
            
            % all T2w_i exprssions should have (1-f_csf) in denominator
            R_dT2_white = (obj.consts.Sm*exp(obj.consts.TE/obj.consts.T2m) ...
                        * obj.Dm * obj.consts.f_white ...
                        * obj.consts.cw_white) ...
                        / (obj.consts.Sw * (obj.consts.f_white + obj.consts.f_grey) * obj.Dw_white);

            val = ((obj.consts.TE * exp(-obj.consts.TE/obj.consts.T2_white)) ...
                / (obj.consts.T2_white^2) ) * R_dT2_white;
        end
        
        % dCm/dT2_csf
        function val = partial_T2csf(obj)
            
            % all T2w_i exprssions should have (1-f_csf) in denominator
            R_dT2_csf = (obj.consts.Sm*exp(obj.consts.TE/obj.consts.T2m) ...
                      * obj.Dm * obj.consts.f_csf ...
                      * obj.consts.cw_csf) ...
                      / (obj.consts.Sw * (obj.consts.f_white + obj.consts.f_grey) * obj.Dw_csf);
            
            val = ((obj.consts.TE * exp(-obj.consts.TE/obj.consts.T2_csf)) ...
                / (obj.consts.T2_csf^2)) * R_dT2_csf;
        end
        
        % dCm/dcw_grey
        function val = partial_cwgrey(obj)
            val = obj.cwi_shared * (obj.consts.f_grey) ...
                / (exp(obj.consts.TE / obj.consts.T2_grey) * obj.Dw_grey);
        end
        
        % dCm/dcw_white
        function val = partial_cwwhite(obj)
            val = obj.cwi_shared * (obj.consts.f_white) ...
                / (exp(obj.consts.TE / obj.consts.T2_white) * obj.Dw_white);
        end
        
        % dCm/dcw_csf
        function val = partial_cwcsf(obj)
            val = obj.cwi_shared * (obj.consts.f_csf) ...
                / (exp(obj.consts.TE / obj.consts.T2_csf) * obj.Dw_csf);
        end
        
        % absolute quantification, Cm
        function val = Cm(obj)
           val = (obj.consts.Sm  .* exp(obj.consts.TE ./ obj.consts.T2m) .* obj.Dm) ...
               ./ (obj.consts.Sw .* (obj.consts.f_white + obj.consts.f_grey)) .* obj.f_sum;
        end
        
        % T1 relaxation weighting by metabolite
        function val = Dm(obj)
            val = (1-exp(-obj.consts.TR./obj.consts.T1m)).^-1;
        end
        
        % T1 relaxation weighting of water in grey matter
        function val = Dw_grey(obj)
            val = (1-exp(-obj.consts.TR./obj.consts.T1_grey)).^-1;
        end
        
        % T1 relaxation weighting of water in white matter
        function val = Dw_white(obj)
            val = (1-exp(-obj.consts.TR./obj.consts.T1_white)).^-1;
        end
        
        % T1 relaxation weighting of water in CSF matter
        function val = Dw_csf(obj)
            val = (1-exp(-obj.consts.TR./obj.consts.T1_csf)).^-1;
        end
        
        % shared term, for water concentration in grey, white CSF
        function val = cwi_shared(obj)
            val = (obj.consts.Sm * exp(obj.consts.TE/obj.consts.T2m) ...
                * obj.Dm) / (obj.consts.Sw * (obj.consts.f_white + obj.consts.f_grey));
        end

        % summation term in Cm
        function val = f_sum(obj)
            val = ((1-obj.consts.f_white-obj.consts.f_grey) .* obj.consts.cw_csf)   ... 
                    ./ exp(obj.consts.TE./obj.consts.T2_csf)   ./ obj.Dw_csf  ...
                + (obj.consts.f_grey  .* obj.consts.cw_grey)  ... 
                    ./ exp(obj.consts.TE./obj.consts.T2_grey)  ./ obj.Dw_grey ...
                + (obj.consts.f_white .* obj.consts.cw_white) ...
                    ./ exp(obj.consts.TE./obj.consts.T2_white) ./ obj.Dw_white ;
        end
        
        % returns set of f_i terms as scalar constants
        function R = shared_fi_terms(obj)
            
            % fi-specific scalar terms
            R.fcsf   = obj.consts.cw_csf / obj.Dw_csf ...
                     / exp(obj.consts.TE / obj.consts.T2_csf);
            R.fgrey  = obj.consts.cw_grey / obj.Dw_grey ...
                     / exp(obj.consts.TE / obj.consts.T2_grey);
            R.fwhite = obj.consts.cw_white / obj.Dw_white ...
                     / exp(obj.consts.TE / obj.consts.T2_white);
                 
            % combined scalar term
            R.fshared = obj.consts.Sm * exp(obj.consts.TE/obj.consts.T2m) ...
                      * obj.Dm/(obj.consts.Sw*(obj.consts.f_white + obj.consts.f_grey)^2);
        end
    end
end