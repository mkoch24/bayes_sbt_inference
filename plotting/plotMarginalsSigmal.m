function [x_mean,x_max_den] = plotMarginalsSigmal(x,x_end,x_min,x_max,...
  k_max,log_scale,plot_true,x_true,varargin)

% returns the mean of the samples of l and sigma^2 corresponding to k_opt
% and the value of l and sigma^2 corresponding to the max. posterior
% density

txt         = varargin{1,1}{1,1};
bw          = varargin{1,1}{1,2};
lw          = varargin{1,1}{1,3};
x0          = varargin{1,1}{1,4};
y0          = varargin{1,1}{1,5};
width       = varargin{1,1}{1,6};
height      = varargin{1,1}{1,7};
col         = varargin{1,1}{1,8};

x_mean      = zeros(k_max,1);
x_max_den   = zeros(k_max,1);
fs          = 12;

figure

for i = 1:k_max

  [den,xi]=ksdensity(x(i,:),'Support',[x_min,x_max],'Bandwidth',bw(i));
  %[den,xi]=ksdensity(x(i,:),'Support',[x_min,x_max]);

  den_max_ind = find(den == max(den));
  x_max_den(i) = xi(den_max_ind(1));  
  
  subplot(1,k_max,i)
  hold on
  plot(xi,den,'LineWidth',lw,'Color',col);
  if plot_true == 1
      plot([x_true(i) x_true(i)],[0 max(den)],'Color',[0.45 0.45 0.45],'LineWidth',1);
  end
  x_mean(i) = mean(x(i,1:x_end));
  
  scatter(x_max_den(i),0,50,'x','MarkerEdgeColor',[0.9 0 0]);

  %figure properties
  axis tight
  if log_scale == 1
    set(gca, 'XScale', 'log')
  end
  xlim([x_min x_max]);
  % xticks([10^-3 10^-1 10^1 10^3]);
  % xticks([10^-2 10^0 10^2]);
  % xticks([10^-3 10^-1 10^1]);
  xlab = sprintf('${%s}_{%i}$',txt,i);
  ylab = sprintf('$p({%s}_{%i}{%s})$',txt,i,['|\bf{m}_{\it{i}},' ...
    '{\it{c}}_{\it{i}},\bf{y}_{\it{i}}']);%'|\bf{m},\bf{\sigma},\bf{y}'
  hXLabel = xlabel(xlab,'Interpreter','latex');
  hYLabel = ylabel(ylab,'Interpreter','latex');
  set( gca                       , ...
    'FontName'   , 'Times New Roman', ...
    'FontSize'   , fs);
  set(hXLabel, ...
    'FontName'   , 'Times New Roman', ...
    'FontSize'   , fs);
  set(hYLabel, ...
    'FontName'   , 'Times New Roman', ...
    'FontSize'   , fs);
  set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'off'      , ...
    'XGrid'       , 'off'      , ...
    'XColor'      , [0 0 0], ...
    'YColor'      , [0 0 0], ...
    'LineWidth'   , 0.5         );  
  hold off

end
set(gcf,'units','centimeters','position',[x0,y0,width,height]);
hold off

end