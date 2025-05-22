function [theta_N_mean] = plotSbt(x,y_obs,y_depth,z_max_pos_prob,...
    k_opt,z_min,z_max,varargin)

% plot sbt figs showing the mean and 1 std-dev envelope

burn_in_samp  = varargin{1,1}{1,1};
n_samp        = varargin{1,1}{1,2};
plot_true     = varargin{1,1}{1,3};
k_true        = varargin{1,1}{1,4};
z_true        = varargin{1,1}{1,5};
theta_true    = varargin{1,1}{1,6};
lw            = varargin{1,1}{1,7};
x_lab         = varargin{1,1}{1,8};
y_lab         = varargin{1,1}{1,9};
x0            = varargin{1,1}{1,10};
y0            = varargin{1,1}{1,11};
width         = varargin{1,1}{1,12};
height        = varargin{1,1}{1,13};

x_min         = min(min(x));
x_max         = max(max(x));
z_true_mod    = [z_min z_true z_max];

figure
hold on
ax1           = gca;
ax1.YDir      = 'reverse';
col           = [0.45,0.45,0.45];
sz            = 20;

scatter(y_obs,y_depth,sz,'x','MarkerEdgeColor',col);
if plot_true == 1
    for i = 1:k_true
        plot([theta_true(i) theta_true(i)],[z_true_mod(i) z_true_mod(i+1)],...
            'Color','k','LineStyle','--', 'Linewidth',lw);
    end
end

theta_N_mean  = mean(x(:,burn_in_samp+1:n_samp),2);
theta_N_std   = std(x(:,burn_in_samp+1:n_samp),0,2);

plot([1.31 1.31],[z_min z_max],'Color',col,'LineStyle',':','LineWidth',lw);
plot([2.05 2.05],[z_min z_max],'Color',col,'LineStyle',':','LineWidth',lw);
plot([2.60 2.60],[z_min z_max],'Color',col,'LineStyle',':','LineWidth',lw);
plot([2.95 2.95],[z_min z_max],'Color',col,'LineStyle',':','LineWidth',lw);
plot([3.60 3.60],[z_min z_max],'Color',col,'LineStyle',':','LineWidth',lw);

ff=0.9;
plot(theta_N_mean,y_depth,'Color',[ff, 0, 0],'LineWidth',lw);
plot(theta_N_mean+theta_N_std,y_depth(1:1:end),'-','Color',[ff, 0.59, 0.59],...
    'LineWidth',lw);
plot(theta_N_mean-theta_N_std,y_depth(1:1:end),'-','Color',[ff, 0.59, 0.59],...
    'LineWidth',lw);

col1 = [0 0 0];
for i = 1:k_opt-1
    plot([x_min(1) x_max(1)],[z_max_pos_prob(i) z_max_pos_prob(i)],...
        'Color',col1,'LineStyle','-','LineWidth',lw);
end

hold off

xlim([x_min(1) x_max(1)]);
xticks([1,2,3,4])
ylim([z_min z_max])
xlab    = sprintf('{%s}',x_lab);
ylab    = sprintf('{%s}',y_lab);
hXLabel = xlabel(xlab);
hYLabel = ylabel(ylab);
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
set(ax1, 'Layer', 'top')


set(gcf,'units','centimeters','Position',[x0,y0,width,height]);
set(gcf, 'PaperPositionMode', 'auto');
pos = ax1.Position;

% create new, empty axes with box but without ticks
ax2 = axes('Position',pos,'box','on','xtick',[],'ytick',[]);
set(ax2,'Color','none');

end