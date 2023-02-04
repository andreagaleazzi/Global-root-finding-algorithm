close all
clear
clc

%% SETUP
x_min = -5.1;
x_max = 5;
n_subdomains = 21;
plot_resolution = 1000;
% i_subdomains = (1:n_subdomains);
grid = linspace(x_min,x_max,n_subdomains)';
% grid = linspace(x_min,x_max,10000)';
plot_grid = linspace(x_min,x_max,plot_resolution)';

tol_fun = 1e-12;
tol_x = 1e-10;
options = optimset('TolFun', tol_fun);
% options = optimset('TolFun', tol_fun, 'Display', 'iter');

f = func(grid);
f_dx = func_dx(grid);

% figure(1)
% hold on
% % plot(grid, f, grid, f_dx);
% % plot(grid, f, 'x', 'MarkerSize', 7);
% plot(grid, f_dx, 'x', 'MarkerSize', 7, 'Color', [0.8500 0.3250 0.0980]);
% plot(grid, f, '--', 'Color', [0 0.4470 0.7410 0.5]);
% plot(grid, f_dx, '-', 'Color', [0.8500 0.3250 0.0980 0.5]);
% % plot(plot_grid, func(plot_grid), ':', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410 0.25]);
% % plot(plot_grid, func_dx(plot_grid), ':', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980 0.25]);
% xline(grid, 'Color', [0,0,0], 'Alpha', 0.1)
% yline(0)
% hold off

boundaries = [
    grid(1:end-1), grid(2:end),...
    f(1:end-1), f(2:end),...
    f_dx(1:end-1), f_dx(2:end)];

% Indexes
i_gridL = 1;
i_gridR = 2;
i_fL = 3;
i_fR = 4;
i_fdxL = 5;
i_fdxR = 6;

%% SOLUTION
domains_with_fdx_roots = boundaries((boundaries(:,i_fdxR) .* boundaries(:,i_fdxL)) <= tol_fun, :);
domains_without_fdx_roots = boundaries((boundaries(:,i_fdxR) .* boundaries(:,i_fdxL)) > tol_fun, :);

domains_with_f_roots_init = boundaries((boundaries(:,i_fR) .* boundaries(:,i_fL)) <= tol_fun, :);
domains_without_f_roots_init = boundaries((boundaries(:,i_fR) .* boundaries(:,i_fL)) > tol_fun, :);

fdx_roots = zeros(size(domains_with_fdx_roots, 1), 2);

for i = 1:size(domains_with_fdx_roots, 1)
    L = domains_with_fdx_roots(i, i_gridL);
    R = domains_with_fdx_roots(i, i_gridR);
    
    try %solving with Bolzano theorem:
        [fdx_roots(i,1), fdx_roots(i,2)] = fzero(@func_dx, [L, R], options);
    catch %if bounds do not change of sign, then numerically relax bounds by tol_fun:
        [fdx_roots(i,1), fdx_roots(i,2)] = fzero(@func_dx, [L-abs(L*tol_x), R+abs(R*tol_x)], options);
    end
end

%%% FUNCTION SEARCH BOUNDARIES
% Defintion of search boundaries on f by mergin fdx solutions with
% initial subdomains with a change of sign on f.
f_boundaries = [fdx_roots(:,1), func(fdx_roots(:,1))];

f_boundaries = uniquetol(f_boundaries,tol_x,'ByRows',true);

% COINCIDENT ROOTS 
f_roots_coincidents = f_boundaries((abs(f_boundaries(:,2)) <= tol_fun), :);
% f_roots_coincidents(:,2) = func(f_roots_coincidents(:,1));

% Merge objective function with initial f domains with roots
f_boundaries = [...
    f_boundaries;...
    domains_with_f_roots_init(:,i_gridL), domains_with_f_roots_init(:,i_fL);...
    domains_with_f_roots_init(:,i_gridR), domains_with_f_roots_init(:,i_fR)...
    ];

% Remove duplicate roots and sort
f_boundaries = uniquetol(f_boundaries,tol_x,'ByRows',true);

% Remove coincident roots
f_boundaries = setdiff(f_boundaries,f_roots_coincidents,'rows');

f_boundaries_with_non_coincident_root = f_boundaries;

domains_with_f_non_coincident_roots = [
    f_boundaries_with_non_coincident_root(1:end-1,1), f_boundaries_with_non_coincident_root(2:end,1),...
    f_boundaries_with_non_coincident_root(1:end-1,2), f_boundaries_with_non_coincident_root(2:end,2),...
    ];

domains_with_f_non_coincident_roots = domains_with_f_non_coincident_roots(...
    (domains_with_f_non_coincident_roots(:,3) .* domains_with_f_non_coincident_roots(:,4)) <= tol_fun, :);

f_roots_non_coincident = zeros(size(domains_with_f_non_coincident_roots, 1), 2);

% domains_with_f_roots

for i = 1:size(domains_with_f_non_coincident_roots, 1)
    L = domains_with_f_non_coincident_roots(i, 1);
    R = domains_with_f_non_coincident_roots(i, 2);
%     x0 = (L+R)/2;
    
    [f_roots_non_coincident(i,1), f_roots_non_coincident(i,2)] = fzero(@func, [L R], options);
end



%% PLOTTING
rgb_f = [0 0.4470 0.7410];
rgb_fdx = [0.8500 0.3250 0.0980];
black = [0 0 0];
white = [1 1 1];

y_min = min(min(boundaries(:,3:end)));
y_max = max(max(boundaries(:,3:end)));

width = 1000;
height = 800;
x0 = (2560 - width)/2;
y0 = (1080 - height)/2;

fig = figure('Name', 'Subplot');
set(gcf,'position',[x0,y0,width,height])

t = tiledlayout(2,2,'TileSpacing','Compact');
% xlabel(t,'x')
% ylabel(t,'y')

%% Analytical plot
% subplot(2,2,1)
nexttile

title('(a) Analytical plot and grid')

hold on

xlim([x_min x_max])
ylim([y_min y_max])



% plot(grid, f, grid, f_dx);
% plot(grid, f, 'x', 'MarkerSize', 7);
% plot(grid, f_dx, 'x', 'MarkerSize', 7, 'Color', [0.8500 0.3250 0.0980]);
% plot(grid, f, '--', 'Color', [0 0.4470 0.7410 0.5]);
% plot(grid, f_dx, '-', 'Color', [0.8500 0.3250 0.0980 0.5]);
plot(plot_grid, func(plot_grid), '-', 'LineWidth', 1, 'Color', rgb_f, 'DisplayName','f(x)');
plot(plot_grid, func_dx(plot_grid), '--', 'LineWidth', 1, 'Color', rgb_fdx, 'DisplayName','fdx(x)');

xline(grid, 'Color', black, 'Alpha', 0.1)
yline(0, 'Color', black, 'Alpha', 0.33)
% legend('$f(x)$', "$f'(x)$", 'Grid', '', 'Interpreter','latex');
legend('f(x)', "f'(x)", 'Grid', '');

hold off


%% Derivative Search
% fig_fdx = figure('Name', 'Derivative search');
% subplot(2,2,2)
nexttile

title('(b) Derivative root finding')

hold on

xlim([x_min x_max])
ylim([y_min y_max])

% Solutions
plot(fdx_roots(:,1), fdx_roots(:,2), 'o', 'Color', black, 'MarkerFaceColor', rgb_fdx)

% plot(grid, f, grid, f_dx);
% plot(grid, f, 'x', 'MarkerSize', 7, 'Color', rgb_f);
% plot(grid, f_dx, '.', 'MarkerSize', 8, 'Color', rgb_fdx);
% plot(grid, f_dx, 'x', 'MarkerSize', 8, 'Color', rgb_fdx);
% plot(grid, f, '--', 'Color', [rgb_f 0.5]);
% plot(grid, f_dx, '-', 'Color', [rgb_fdx 0.5]);
% plot(plot_grid, func(plot_grid), ':', 'LineWidth', 1.5, 'Color', [rgb_f 0.25]);
% plot(plot_grid, func_dx(plot_grid), ':', 'LineWidth', 1.5, 'Color', [rgb_fdx 0.25]);


for i = 1:size(domains_without_fdx_roots(:,1), 1)
    plot(domains_without_fdx_roots(i, i_gridL:i_gridR), domains_without_fdx_roots(i, i_fdxL:i_fdxR), ':', 'LineWidth', 1.5, 'Color', [black 0.25])
end
plot(domains_without_fdx_roots(:, i_gridL), domains_without_fdx_roots(:, i_fdxL), 'x', 'MarkerSize', 8, 'Color', [black 0.25]);
plot(domains_without_fdx_roots(:, i_gridR), domains_without_fdx_roots(:, i_fdxR), 'x', 'MarkerSize', 8, 'Color', [black 0.25]);

for i = 1:size(domains_with_fdx_roots(:,1), 1)
    plot(domains_with_fdx_roots(i, i_gridL:i_gridR), domains_with_fdx_roots(i, i_fdxL:i_fdxR), 'LineWidth', 1.5, 'Color', [rgb_fdx 0.8])
end
plot(domains_with_fdx_roots(:, i_gridL), domains_with_fdx_roots(:, i_fdxL), 'x', 'MarkerSize', 8, 'Color', rgb_fdx);
plot(domains_with_fdx_roots(:, i_gridR), domains_with_fdx_roots(:, i_fdxR), 'x', 'MarkerSize', 8, 'Color', rgb_fdx);

% plot(plot_grid, func_dx(plot_grid), '-.', 'LineWidth', 1, 'Color', [rgb_fdx, 0.25], 'DisplayName','fdx(x)');

% plot(f_roots(:,1), f_roots(:,2), 'o', 'Color', rgb_f)
% plot(f_roots_coincidents(:,1), f_roots_coincidents(:,2), 'o', 'Color', [rgb_f 0.5], 'MarkerFaceColor', rgb_f)


xline(grid, 'Color', black, 'Alpha', 0.1)
yline(0, 'Color', black, 'Alpha', 0.33)

% legend('', '', '', '', '', '','f(x) root', "f'(x) root", "Coincident root", '');
% legend('color', 'w')
% legend('', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '')
% legend("f'(x) root", 'test', '')
legend("f'(x) root", '','','','','','','','','','','','','','','','','','','','','Node value')

hold off

%% Objective Search
% fig_f = figure('Name', 'Objective search');
% subplot(2,2,3)
nexttile
hold on

title('(c) Function non-coincident root finding')

set(gca,'Color', [black 0.3])
xlim([x_min x_max])
ylim([y_min y_max])


for i = 1:1:size(domains_with_f_non_coincident_roots, 1)
% for i = 2
 	rectangle('Position',[domains_with_f_non_coincident_roots(i,1), y_min, domains_with_f_non_coincident_roots(i,2)-domains_with_f_non_coincident_roots(i,1), y_max-y_min],...
         'FaceColor',  white, 'EdgeColor', [0 0 0 0], 'LineWidth', 1)
end
yline(0, 'Color', black, 'Alpha', 0.33)

for i = 1:1:size(domains_with_f_non_coincident_roots(:,1), 1)
    plot(domains_with_f_non_coincident_roots(i, 1:2), domains_with_f_non_coincident_roots(i, i_fL:i_fR), 'LineWidth', 1.5, 'Color', [rgb_f 0.8])
end
plot(domains_with_f_non_coincident_roots(:, i_gridL), domains_with_f_non_coincident_roots(:, i_fL), 'x', 'MarkerSize', 8, 'Color', rgb_f);
plot(domains_with_f_non_coincident_roots(:, i_gridR), domains_with_f_non_coincident_roots(:, i_fR), 'x', 'MarkerSize', 8, 'Color', rgb_f);

% Solutions
plot(f_roots_non_coincident(:,1), f_roots_non_coincident(:,2), 'o', 'Color', black, 'MarkerFaceColor', rgb_f)

plot(fdx_roots(:,1), fdx_roots(:,2), 'o', 'Color', black, 'MarkerFaceColor', rgb_fdx)

plot(f_roots_coincidents(:,1), f_roots_coincidents(:,2), 'square', 'MarkerSize', 8, 'Color', black, 'MarkerFaceColor', black)

plot(plot_grid, func(plot_grid), ':', 'LineWidth', 1.5, 'Color', [black 0.25]);
xline(f_boundaries(:,1), 'Color', black, 'Alpha', 0.1)
yline(y_min, 'Color', black)

legend('', '', '', '', '', 'f(x) root', '', 'Coincident root', 'f(x)', '');
legend('color', 'w')

hold off


%% Solutions plot
% subplot(2,2,4)
nexttile

hold on

title('(d) Function roots collection')
xlim([x_min x_max])
ylim([y_min y_max])



% plot(grid, f, grid, f_dx);
% plot(grid, f, 'x', 'MarkerSize', 7);
% plot(grid, f_dx, 'x', 'MarkerSize', 7, 'Color', [0.8500 0.3250 0.0980]);
% plot(grid, f, '--', 'Color', [0 0.4470 0.7410 0.5]);
% plot(grid, f_dx, '-', 'Color', [0.8500 0.3250 0.0980 0.5]);
plot(plot_grid, func(plot_grid), '-', 'LineWidth', 1, 'Color', rgb_f, 'DisplayName','f(x)');
plot(plot_grid, func_dx(plot_grid), '--', 'LineWidth', 1, 'Color', rgb_fdx, 'DisplayName','fdx(x)');

% xline(grid, 'Color', black, 'Alpha', 0.1)
yline(0, 'Color', black, 'Alpha', 0.33)


plot(f_roots_non_coincident(:,1), f_roots_non_coincident(:,2), 'o', 'Color', black, 'MarkerFaceColor', rgb_f)
plot(f_roots_coincidents(:,1), f_roots_coincidents(:,2), 'o', 'Color', black, 'MarkerFaceColor', rgb_f)

% legend('$f(x)$', "$f'(x)$", '', 'Interpreter','latex');

hold off

%% TWEAKS

% set(findall(gca,'-property','FontName'),'FontName',Times New Roman)
% set(findall(gcf,'-property','FontSize'),'FontSize',10)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

%% EXPORT FIGURE

export_fig(fig, 'fig', 'subplot')