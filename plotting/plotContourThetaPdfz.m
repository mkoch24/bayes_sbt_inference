function plotContourThetaPdfz(nor,x,y,nx,ny,y_depth,varargin)

% plot the contour of theta and the marginal posterior pdf of interface
% depths considering all possible models sampled in the algorithm

%% ========================================================================
% plot of the contour of theta

x_min = min(x);
x_max = max(x);

figure
subplot(1,3,[1,2])
hold on
histogram2(x,y,[nx ny],'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization',nor); colorbar("westoutside");
ax1 = gca;
ax1.YDir = 'reverse';
col = [0.45,0.45,0.45];
sz = 20;

if varargin{1,1}{1,1} == 1
    % add confidence intervals

    d             = varargin{1,1}{1,2};
    x_mat         = varargin{1,1}{1,3};
    burn_in_samp  = varargin{1,1}{1,4};
    n_samp        = varargin{1,1}{1,5};
    k_opt         = varargin{1,1}{1,6};
    z_choice      = varargin{1,1}{1,7};
    y_obs         = varargin{1,1}{1,8};
    lw            = varargin{1,1}{1,9};
    x_lab         = varargin{1,1}{1,10};
    y_lab         = varargin{1,1}{1,11};
    z_min         = varargin{1,1}{1,12};
    z_max         = varargin{1,1}{1,13};
    x0            = varargin{1,1}{1,14};
    y0            = varargin{1,1}{1,15};
    width         = varargin{1,1}{1,16};
    height        = varargin{1,1}{1,17};

    
    scatter(y_obs,y_depth,sz,'x','MarkerEdgeColor',col);

    theta_N_mean = mean(x_mat(:,burn_in_samp+1:n_samp),2);
    theta_N_sort = zeros(d,n_samp-burn_in_samp);
    for i = 1:d
        theta_N_sort(i,:) = sort(x_mat(i,burn_in_samp+1:n_samp),2);
    end
    theta_N_5  = theta_N_sort(:,round((n_samp-burn_in_samp)*0.025));
    theta_N_50 = theta_N_sort(:,round((n_samp-burn_in_samp)*0.5));
    theta_N_95 = theta_N_sort(:,round((n_samp-burn_in_samp)*0.975));

    %%{
    ff=0.9;
    plot(theta_N_mean,y_depth,'Color',[ff, 0, 0],'LineWidth',lw);
    plot(theta_N_5(1:1:end),y_depth(1:1:end),'-','Color',[ff, 0.59, 0.59],...
        'LineWidth',lw);
    plot(theta_N_50(1:1:end),y_depth(1:1:end),'-','Color',[ff/2, 0, 0],...
        'LineWidth',lw);
    plot(theta_N_95(1:1:end),y_depth(1:1:end),'-','Color',[ff, 0.59, 0.59],...
        'LineWidth',lw);

    %%{
    col1 = [0 0 0];
    for i = 1:k_opt-1
        plot([x_min(1) x_max(1)],[z_choice(i) z_choice(i)],'Color',col1,...
            'LineStyle','-','LineWidth',lw/2);
    end
    %}
end
hold off

xlim([x_min(1) x_max(1)]);
xticks([1,2,3,4])
ylim([z_min z_max])
xlab = sprintf('{%s}',x_lab);
ylab = sprintf('{%s}',y_lab);
hXLabel = xlabel(xlab);
hYLabel = ylabel(ylab);
% colorbar
colormap sky

fs = 12;
set( ax1                       , ...
    'FontName'   , 'Times New Roman', ...
    'FontSize'   , fs);
set(hXLabel, ...
    'FontName'   , 'Times New Roman', ...
    'FontSize'   , fs);
set(hYLabel, ...
    'FontName'   , 'Times New Roman', ...
    'FontSize'   , fs);

set(ax1, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.01 .01] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [0 0 0], ...
    'YColor'      , [0 0 0], ...
    'LineWidth'   , 0.25         );
set(ax1, 'Layer', 'top')

%% ========================================================================
% plot of the posterior marginal pdf of z

theta_min       = varargin{1,1}{1,18};
theta_max       = varargin{1,1}{1,19};
int_depth_col   = varargin{1,1}{1,20};
col_dens        = varargin{1,1}{1,21};
delta_int       = varargin{1,1}{1,22};
r               = varargin{1,1}{1,23};
plot_true       = varargin{1,1}{1,24};
k_true          = varargin{1,1}{1,25};
z_true          = varargin{1,1}{1,26};
theta_true      = varargin{1,1}{1,27};


subplot(1,3,3); 

% Get the position of this subplot
pos = get(gca, 'Position'); 
delete(gca);  % Delete default axes to replace with custom axes

col     = [0.8,0.8,0.8];
sz      = 20;

% Primary axis
ax1 = axes('Position', pos);
hold(ax1, 'on');
scatter(ax1, y_obs, y_depth, sz, 'x', 'MarkerEdgeColor', col);
z_true_mod = [z_min z_true z_max];
if plot_true == 1
    for i = 1:k_true
        plot([theta_true(i) theta_true(i)],[z_true_mod(i) z_true_mod(i+1)],...
            'Color','k','LineStyle','--', 'Linewidth',lw);
    end
end
ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';
ax1.YDir = 'reverse';
ax1.YLim = [z_min z_max];
ax1.XLim = [theta_min(1) theta_max(1)];
ax1.YTickLabel = [];
bw = 0.05;
xticks(ax1, [1, 2, 3, 4]);
x_lab1 = '\it I_{c}';
xlab1 = sprintf('{%s}', x_lab1);
hXLabel1 = xlabel(ax1, xlab1);

% Secondary axis
ax2 = axes('Position', pos);
[den_z, xi_z] = ksdensity(int_depth_col, y_depth, 'Bandwidth', bw);
plot(ax2, den_z, xi_z, 'Color', col_dens, 'LineWidth', lw);
ax2.YDir = 'reverse';
ax2.XAxisLocation = 'top';
ax2.Color = 'none';
set(ax2, 'YTick', []);
ax2.YLim = [z_min + r * delta_int, z_max - r * delta_int];
hXLabel2 = xlabel(ax2, '$p(z|\bf{y})$', 'Interpreter', 'latex');

% Formatting
fs = 12; % Example font size
set([ax1, ax2], 'FontName', 'Times New Roman', 'FontSize', fs);
set(hXLabel1, 'FontName', 'Times New Roman', 'FontSize', fs);
set(hXLabel2, 'FontName', 'Times New Roman', 'FontSize', fs);
set(ax2, 'Layer', 'top');
hold(ax1, 'off');
hold(ax2, 'off');


set(gcf,'units','centimeters','Position',[x0,y0,width,height]);
set(gcf, 'PaperPositionMode', 'auto');

end