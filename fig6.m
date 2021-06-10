% This file generates Figure 6 of J. Fluid Mech. (2018), vol. 854, pp.
% 34-55

% clear workspace
clear; close all

% add support files
addpath('external_files')

% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,sqrt(2),7/4]*monitor_size(4)/2)

% margins for figure
xl = 0.10; % left margin
xr = 0.03; % right margin
xg = 0.14; % gap between plots
xa = (1 - xl - xr - sum(xg))/ (length(xg) + 1); % width of each plot
yt = 0.01; % top margin
yb = 0.05; % bottom margin
yg = [0.055, 0.055];
ya = (1 - yt - yb - sum(yg)) / (length(yg) + 1); % height of each plot

% x coordinates for subfigures
px(1) = xl;
px(2) = px(1) + xa + xg;
px([3:4,5]) = px([1,1,2]);

% y coordinate for subfigures
py(4) = yb;
py(3) = py(4) + ya + yg(2);
py(1) = py(3) + ya + yg(1);
py([5;2]) = py([3,1]);

ffsize = 13; % font size of figure
flsize = get(0,'DefaultLineLinewidth')*1; % linewidth of figure
fmsize = get(0,'DefaultLineMarkerSize')*1.4; % marker size of figure

% set syles for figures
mu0 = [0.41, 0.56, 0.71];
stylename = {'k-','k--','k-.','k:'};
stylename2 = {'ko','kd','ksq'};
contour_vec = (1:log10(10^(1/4)):20);

% files to load in
filename = {'LQG_041.mat','LQG_056.mat','LQG_071.mat'};
%% a - c

% loop through plots
for i = 1:length(filename)
    
    % load in data
    load(['fig6/',filename{i}])
    
    % create subplot
    subplot('position', [px(i) py(i) xa ya])
    
    x1 = - sqrt( 2*mu0(i) / 0.01);
    x2 = -x1;
    
    % draw background
    fill([x1 x2 x2 x1 x1], [x1 x1 x2 x2 x1],[1 1 1] * 0.95,'edgecolor',[0.9 0.9 0.9])
    hold on 
    plot([x1 x2 x2 x1 x1], [x1 x1 x2 x2 x1],'k','LineWidth',flsize,'LineStyle',':')

    % find optimals
    [mag,a] = min(gamma_2_IOC.');
    [~,b] = min(mag);
    
    xs_vec = (x_vec);
    xa_vec = (x_vec);
    
    c = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10];
    contour(xs_vec,xa_vec,log10(real(gamma_2_IOC)).',contour_vec,'linewidth',flsize);
    caxis(log10([c(1) c(length(c))]));
    
    plot (xa_vec(b),xs_vec(a(b)),stylename2{i},'markersize',fmsize,'linewidth',flsize);
    set(gca,'FontSize',ffsize)
    
    
    
    xlabel('$x_s$','interpreter','latex','FontSize',ffsize)
    set(gca,'XTick',-15:5:15);
    
    ylabel('$x_a$','interpreter','latex','FontSize',ffsize)
    
    
    if i == 1
        text(-20.93,15,'$(a)$','interpreter','latex','FontSize',ffsize)
    elseif i == 2
        text(-20.93,15,'$(b)$','interpreter','latex','FontSize',ffsize)
    elseif i == 3
        text(-20.93,15,'$(c)$','interpreter','latex','FontSize',ffsize)
    end
    
    set(gca,'TickLabelInterpreter','Latex')
    caxis([1,11])
    %axis equal
end
colormap(brewermap([],'YlGnBu')), hold on,
colormap(flipud(colormap))

%%






subplot('position', [px(4) py(4) (xa*2 + xg) ya]), hold on
set(gca,'FontSize',ffsize)

load fig6/LQG_mu0range.mat
%%

% define upstream and downstream branch
x12branch = real(sqrt( 2*mu0vec / 0.01));
% cut off last few points of branch  to avoid overlap with box
cutoff = 21;
% plot branches
fill([x12branch(1:end-cutoff), fliplr(-x12branch(1:end-cutoff))],[mu0vec(1:end-cutoff),fliplr(mu0vec(1:end-cutoff))],[0.95 0.95 0.95],'EdgeColor',[0.95,0.95,0.95])
plot([x12branch, NaN, fliplr(-x12branch)],[mu0vec,NaN,fliplr(mu0vec)],stylename{4})

plot(x_s_vec(1:1:end),mu0vec(1:1:end),'k-.','linewidth',flsize)
for i = 1:length(mu0)
    plot(x_s_vec(mu0vec == mu0(i)),mu0vec(mu0vec == mu0(i)),stylename2{i},'markersize',fmsize,'linewidth',flsize);
end
%%
ylim([-0.01,0.90])
set(gca,'YTick',[0:0.1:0.9],'XTick',-15:2.5:15)
%%
plot(x_a_vec(1:1:end),mu0vec(1:1:end),'k:','linewidth',flsize)%,'o','color','r','markersize',4)
for i = 1:length(mu0)
    plot(x_a_vec(mu0vec == mu0(i)),mu0vec(mu0vec == mu0(i)),stylename2{i},'markersize',fmsize,'linewidth',flsize);
end%%
ylim([-0.01,0.9])
ylabel('$$ \mu_0 $$','Interpreter','latex','FontSize',ffsize),
xlabel(' $$ x $$','Interpreter','latex','FontSize',ffsize)%xlabel('$$ LQG: \{x_{s-opt}(-\cdot), x_{a-opt}(:)\}$$ and $$LQE/LQR: \{x_{s-opt}(-), x_{a-opt}(--)\} $$','Interpreter','latex')
set(gca,'YTick',[0:0.1:0.9],'XTick',-15:2.5:15)
box on
%%
%%

load fig4/LQE_LQR_mu0range.mat
%%
plot(x_s_vec(1:1:end),mu0vec(1:1:end),stylename{1},'linewidth',flsize)
plot(x_a_vec(1:1:end),mu0vec(1:1:end),stylename{2},'linewidth',flsize)

text(-17.7,0.9,'$(d)$','Interpreter','latex','FontSize',ffsize)
set(gca,'TickLabelInterpreter','Latex')






load fig6/LQG_mu0range.mat
%%
subplot('position', [px(5) py(5) xa ya])


semilogy(mu0vec,sqrt(J_vec),'k'), hold on
set(gca,'FontSize',ffsize)
%
for i = 1:length(mu0)
    plot(mu0vec(mu0vec == mu0(i)),sqrt(J_vec(mu0vec == mu0(i))),stylename2{i},'markersize',fmsize,'linewidth',flsize);
end
%%
axis([-0.01,0.9,10^0,1.2*10^3])
xlabel('$$ \mu_0 $$','Interpreter','latex','FontSize',ffsize), ylabel('$$ \gamma_{IO} (x_{a-opt},x_{s-opt}) $$','Interpreter','latex','FontSize',ffsize)

text(-0.19,1.2e3,'$(e)$','Interpreter','latex','FontSize',ffsize)
set(gca,'TickLabelInterpreter','Latex')

% print('-depsc2','fig6')

