function plotChainHistk(burn_in_samp,n_samp,k,k_prior_y,k_min,k_max,varargin)

plot_true   = varargin{1,1}{1,1};
k_true      = varargin{1,1}{1,2};
lw          = varargin{1,1}{1,3};
x0          = varargin{1,1}{1,4};
y0          = varargin{1,1}{1,5};
width       = varargin{1,1}{1,6};
height      = varargin{1,1}{1,7};
col_hist    = varargin{1,1}{1,8};

figure
% plot Markov chain for k
subplot(1,2,1)
plot(k);
ax1         = gca;
ax1.YLim    = [2 k_max];
ax1.YTick   = [0:2:k_max];
hXLabel1    = xlabel('Sample number');
hYLabel1    = ylabel('\it k');

% plot histogram of k
subplot(1,2,2); hold on;
histogram(k(1,burn_in_samp+1:n_samp),'FaceColor',col_hist,...
    'FaceAlpha',1,'EdgeColor','none','Normalization','pdf');
plot([k_min,k_max],[k_prior_y k_prior_y],'LineStyle','--',...
    'LineWidth',lw*2,'Color',[0.45 0.45 0.45]);
if plot_true == 1
    plot([k_true k_true],[0 1],'Color',[0.45 0.45 0.45],...
        'LineWidth',lw*2);
end

ax2 = gca;
ax2.XLim = [2 k_max];
ax2.XTick = [2:6:k_max];
hXLabel2 = xlabel('\it k');
hYLabel2 = ylabel('$p(k|\bf{y})$','Interpreter','latex');
hold off

% figure properties
ax      = [ax1 ax2];
hXLabel = [hXLabel1 hXLabel2];
hYLabel = [hYLabel1 hYLabel2];
fs      = 12;
set( ax                       , ...
    'FontName'   , 'Times New Roman', ...
    'FontSize'   , fs);
set(hXLabel, ...
    'FontName'   , 'Times New Roman', ...
    'FontSize'   , fs);
set(hYLabel, ...
    'FontName'   , 'Times New Roman', ...
    'FontSize'   , fs);
set(ax, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.01 .01] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [0 0 0], ...
    'YColor'      , [0 0 0], ...
    'LineWidth'   , 0.5         );


set(gcf,'units','centimeters','Position',[x0,y0,width,height]);
set(gcf, 'PaperPositionMode', 'auto');

end