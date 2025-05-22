function [z_max_pos_prob,k_max_pos_prob] =postProcessMaxDensz(post_prob,...
    z_ind,z_ind_temp,k_index,k_opt,varargin)

% find and plot the pdf of interface depths corresponding to the sample
% with max conditional posterior density

plot_true   = varargin{1,1}{1,1};
z_true      = varargin{1,1}{1,2};
txt         = varargin{1,1}{1,3};
bw          = varargin{1,1}{1,4};
lw          = varargin{1,1}{1,5};
x0          = varargin{1,1}{1,6};
y0          = varargin{1,1}{1,7};
width       = varargin{1,1}{1,8};
height      = varargin{1,1}{1,9};
col_dens    = varargin{1,1}{1,10};

post_prob_max_index_master = ...
    find(post_prob(z_ind_temp{k_index == k_opt,1},1) == ...
    max(post_prob(z_ind_temp{k_index == k_opt,1},1)));
z_max_pos_prob = ...
    z_ind{k_index == k_opt,1}(1:k_opt-1,post_prob_max_index_master(1));
k_max_pos_prob = length(z_max_pos_prob)+1;
plotMaxDensz(z_ind{k_index == k_opt,1},z_true,k_opt-1,txt,plot_true,lw,...
    col_dens,bw,z_max_pos_prob,x0,y0,width,height);

end


%% ========================================================================
function [x_mean,x_max_den] = plotMaxDensz(x,x_true,kmax,txt,plot_true,...
    lw,col,bw,z_max_pos_prob,x0,y0,width,height)

% plot the posterior density functions of interface depths corresponding to
% the sample with max conditional posterior density  

x_mean      = zeros(kmax,1);
x_max_den   = zeros(kmax,1);

% plot 1-D marginals of z
fs = 12;
figure
for i = 1:kmax

    [den,xi]=ksdensity(x(i,:),'Bandwidth',bw);
    %[den,xi]=ksdensity(x(i,:));
    den_max_ind = find(den == max(den));
    x_max_den(i) = xi(den_max_ind(1));
    subplot(1,kmax,i)
    hold on
    plot(xi,den,'LineWidth',lw,'Color',col);
    if plot_true == 1
        plot([x_true(i) x_true(i)],[0 max(den)],'Color',...
            [0.45 0.45 0.45],'LineWidth',1);
    end

    % plot interface depth corresponding to max conditional posterior prob
    scatter(z_max_pos_prob(i),0,50,'x','MarkerEdgeColor',[0.9 0 0]);

    % figure properties
    axis tight
    xlab = sprintf('${%s}_{%i}$',txt,i);
    %ylab = sprintf('$p({%s}_{%i}{%s})$',txt,i,'|k={%i},\bf{y}',k_opt);
    ylab = sprintf('$p({%s}_{%i}{|k=%i,\\bf{y}})$', txt, i, kmax+1);
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
        'TickLength'  , [.01 .01] , ...
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