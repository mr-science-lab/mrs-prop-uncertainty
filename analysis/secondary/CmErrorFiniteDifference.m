% ------------------------------------------------------------
% Script:  CmErrorFiniteDifference.m
% Summary: approximates the partial derivative of parameters used in
%          absolute quantification of MR spectroscopy datasets using the
%          finite difference approximation. 
% 
% Usage:   >> CmErrorFiniteDifference
%
% Inputs:  None
%
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
%
% Columbia University, 2023
% ------------------------------------------------------------
function CmErrorFiniteDifference(config)
    % verify Partial Derivatives, Analytical Solutions
    % compare approximation to analytical derivation
    close all; clc;

    % set parameters
    n        = config.CmErrorFiniteDifference.n;
    params   = config.CmErrorFiniteDifference.parameters;
    partials = config.CmErrorFiniteDifference.partials;
    xrange   = { [.05 .5],[0 2],[0.05 .5],          ...
                 [0,0.5],[0,0.5],[0 2],[0 2],[0 2], ...
                 [0,0.5], [0,0.5], [0,0.5],         ...
                 [0 100], [0 100], [0 100],         ...
                 [0 0.5], [0.2 2] };

    % run approximation for each parameter
    for i = 1:length(params)

        % set-up
        param = params{i}; % current parameter
        ds    = CmUtils.load_dataset('CmErrorParams.json');  % reset constants

        % set param based on selected range
        ds.constants.(param) = linspace(xrange{i}(1),xrange{i}(2),n);
        range_vals = ds.constants.(param);
        cmobj      = CmPartials(ds.constants);

        % absolute quantification
        Cm = cmobj.Cm;

        % analytical solution
        analyt_soln = zeros(1,n);
        for j = 1:n
            ds.constants.(param) = range_vals(j);
            cmobj = CmPartials(ds.constants);
            analyt_soln(j) = cmobj.(partials{i});
        end

        % finite difference approx
        h     = (xrange{i}(2)-xrange{i}(1))/n;
        fdiff = diff(Cm)./h;

        % plot on single figure
        fig = figure;
        fig.Position = [100 100 550 450];
        set(fig, 'Color', [1 1 1]);
        yyaxis left;
        plot( range_vals,Cm*1000, ...
              'color',[.6 .6 .6], ...
              'LineStyle','-',    ...
              'LineWidth',2);
        hold on;
        grid on;
        xlabel(CmUtils.format_param_as_string(param));
        ylabel('C_m (mM)');
        xlim([xrange{i}(1),xrange{i}(2)]);

        yyaxis right;
        plot( range_vals,analyt_soln, ...
              'color','#0076BA',      ...
              'LineStyle','-',        ...
              'LineWidth',3 );
        hold on;
        plot(range_vals(1:(n-1)),fdiff,'ko','LineWidth',1.5);
        ylabel('\partial{C_m}/\partial{x_i}');
        if strcmp(param,'cw_csf') || strcmp(param,'cw_grey') || strcmp(param,'cw_white')
            ylim([0 1e-4]);
        end

        % format and style
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';
        ax.FontName = 'Calibri';
        ax.FontSize = 18;
        leg = legend('C_m',strcat('\partial{C_m}/\partial{',...
            CmUtils.format_param_as_string(param),'}'),'\Delta[C_m]',...
            'Location','best');
        leg.FontSize = 16;
        leg.FontName = 'Calibri';
        title( CmUtils.format_param_as_string(param),...
               'FontSize',22,        ...
               'FontName','Calibri', ...
               'FontWeight','normal' );

        % save image
        exportgraphics(fig,strcat(config.paths.res_dir,'finite-difference-approx',param,'.png'),'Resolution',2000);
    end
end
