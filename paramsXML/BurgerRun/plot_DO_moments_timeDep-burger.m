% Load Data
[x] = mmread('MeanSolnFinal.mm');
[V] = mmread('V_BASE.mm');
[Y] = mmread('YTrans_COEFF.mm');
xx = linspace(0,1,128);
lw = 1.5;

load Exact_data.mat;
exact_sol=val_0_8;
fig = figure('Units', 'centimeters', 'Position', [2, 2, 24, 18]); % Replaces w and h from tight_subplot
tlo = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Pre-calculate data for the loop
vy = V * Y';
moments_data = {exact_sol.mean, exact_sol.variance, exact_sol.m3, exact_sol.kurtosis_excess};
titles = {'mean', 'variance', '3rd moment', 'excess 4th moment'};
y_labels = {'$\mathbf{E\mathit u}$', '$\mathbf{E\mathit u^2}$', '$\mathbf{E\mathit u^3}$', '$\mathbf{E\mathit u^4-3}(\mathbf{E} \mathbf{\mathit u^2)^2}$'};

% Loop to create each tile
for i = 1:4
    ax = nexttile;
    hold(ax, 'on');
    
    % Plot Exact Solution
    plot(xx', moments_data{i}, 'linewidth', lw, 'DisplayName', 'Exact');
    
    % Calculate and Plot DO Solution
    if i == 1
        plot(xx', x, 'linewidth', lw, 'DisplayName', 'DO');
        legend('Interpreter', 'latex', 'Location', 'best', 'FontSize',12);
    else
        if i == 4
            % Reference mom2 from i=2 step
            mom = sum(vy.^i, 2) / size(Y, 1) - 3 * mom2.^2;
        else
            mom = sum(vy.^i, 2) / size(Y, 1);
            if i == 2, mom2 = mom; end
        end
        plot(xx', mom, 'linewidth', lw, 'DisplayName', 'DO');
    end
    
    % Formatting
    title(titles{i},'FontSize', 14);
    xlabel('$\mathbf{\mathit x}$', 'interpreter', 'latex', 'FontSize', 14);
    ylabel(y_labels{i}, 'interpreter', 'latex', 'FontSize', 14);
    
    if i == 4
        yticks((-2:1:1)*10^-5);
    end
    
    grid(ax, 'on');
end
% Apply to all axes in the current figure
set(findall(fig, 'type', 'axes'), 'Box', 'on', 'XGrid', 'off', 'YGrid', 'off');

% Optional: Set overall figure properties
%set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
