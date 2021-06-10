% This file generates Figure 4 of J. Fluid Mech. (2018), vol. 854, pp.
% 34-55

% clear workspace
clear; close all


% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,sqrt(2),7/4]*monitor_size(4)/2)

% margins for figure
xl = 0.08; % left margin
xr = 0.015; % right margin
xg = 0.10; % gap between plots
xa = (1 - xl - xr - sum(xg))/ (length(xg) + 1); % width of each plot
yt = 0.02; % top margin
yb = 0.06; % bottom margin
yg = [0.08, 0.08];
ya = (1 - yt - yb - sum(yg)) / (length(yg) + 1); % height of each plot

% x coordinates for subfigures
px(1) = xl;
px(2) = px(1) + xa + xg;
px(4:5) = px(1:2);
px(3) = xl;

% y coordinates for subfigures
py(5) = yb;
py(3) = py(5) + ya + yg(2);
py(1) = py(3) + ya + yg(1);
py([4;2]) = py([5,1]);

ffsize = 13; % font size of figure
flsize = get(0,'DefaultLineLinewidth')*1; % linewidth of figure
fmsize = get(0,'DefaultLineMarkerSize')*1.4; % marker size of figure

% set linestyle for figures
stylename = {'k-','k--','k-.','k:'};
stylename2 = {'ko','kd','ksq'};

% files to load in
filename = {'LQE_LQR_041.mat','LQE_LQR_056.mat','LQE_LQR_071.mat'};

%% gamma OE figure (a)
% first subfigure
subplot('position', [px(1) py(1) xa ya])

for jk = 1:3
% load file
load(['fig4/',filename{jk}])
% plot performance parameter
semilogy(x_vec,(gamma_2_OE),stylename{jk},'linewidth',flsize),  hold on
% find minimum
[c,d] = min((gamma_2_OE));
% plot minimum
plot(x_vec(d),(c),stylename2{jk},'markersize',fmsize)
end
% axis properties
axis([-15,15,1,10^5])
xlabel('$x_s$','interpreter','latex','FontSize',ffsize)
ylabel('$\gamma_{OE}$','interpreter','latex','FontSize',ffsize)
text(-20.75,1e5,'$(a)$','interpreter','latex','FontSize',ffsize)
set(gca,'YTick',[1,10,100,1e3,1e4,1e5,1e6],'XTick',-15:5:15,'FontSize',ffsize)
set(gca,'TickLabelInterpreter','Latex')

%% gamma FI figure (b)
% call second subfigure
subplot('position', [px(2) py(2) xa ya])

for jk = 1:3
% repeat previous loop
load(['fig4/',filename{jk}])
semilogy(x_vec,(gamma_2_FI),stylename{jk},'linewidth',flsize),  hold on
[c,d] = min((gamma_2_FI));
plot(x_vec(d),(c),stylename2{jk},'markersize',fmsize)
end

% axis properties
axis([-15,15,1,10^5])
xlabel('$x_a$','interpreter','latex','FontSize',ffsize)
ylabel('$\gamma_{FI}$','interpreter','latex','FontSize',ffsize)
text(-20.75,1e5,'$(b)$','interpreter','latex','FontSize',ffsize')
set(gca,'YTick',[1,10,100,1e3,1e4,1e5,1e6],'XTick',-15:5:15,'FontSize',ffsize)
set(gca,'TickLabelInterpreter','Latex')


%% x_s x_a opt

% third subfigure
subplot('position', [px(3) py(3) (2*xa+xg) ya]), hold on

% load data
load fig4/LQE_LQR_mu0range.mat

% define upstream and downstream branch
x12branch = real(sqrt( 2*mu0vec / 0.01));
% cut off last few points of branch  to avoid overlap with box
cutoff = 21;
% plot branches
fill([x12branch(1:end-cutoff), fliplr(-x12branch(1:end-cutoff))],[mu0vec(1:end-cutoff),fliplr(mu0vec(1:end-cutoff))],[0.95 0.95 0.95],'EdgeColor',[0.95,0.95,0.95])
plot([x12branch, NaN, fliplr(-x12branch)],[mu0vec,NaN,fliplr(mu0vec)],stylename{4})

%% plot x_s

% plot values
plot(x_s_vec(1:1:end),mu0vec(1:1:end),stylename{1},'linewidth',flsize)
% highlight values from paper
mu0 = [0.41, 0.56, 0.71];
stylename2 = {'ko','kd','ksq'};
for jk = 1:length(mu0)
plot(x_s_vec(mu0vec == mu0(jk)),mu0vec(mu0vec == mu0(jk)),stylename2{jk},'markersize',fmsize,'linewidth',flsize);
end

% axis labels
ylabel('$$ \mu_0 $$','Interpreter','latex','FontSize',ffsize), 
xlabel('$$ x_{s-opt} $$','Interpreter','latex','FontSize',ffsize)
ylim([-0.01,0.90])
set(gca,'YTick',[0:0.1:0.9],'XTick',-15:2.5:15,'FontSize',ffsize)
set(gca,'TickLabelInterpreter','Latex')
%


%% plot x_a

% plot values
plot(x_a_vec(1:1:end),mu0vec(1:1:end),stylename{2},'linewidth',flsize)

% highlight values from paper
for jk = 1:length(mu0)
plot(x_a_vec(mu0vec == mu0(jk)),mu0vec(mu0vec == mu0(jk)),stylename2{jk},'markersize',fmsize,'linewidth',flsize);
end

% axis labels
ylim([-0.02,0.9])
ylabel('$$ \mu_0 $$','Interpreter','latex','FontSize',ffsize), 
xlabel('$$ x $$','Interpreter','latex','FontSize',ffsize)%xlabel('$$ \{x_{s-opt}(-), x_{a-opt}(--)\} $$','Interpreter','latex')
set(gca,'YTick',[0:0.1:0.9],'XTick',-15:2.5:15,'FontSize',ffsize)
text(-17.4917,0.9,'$(c)$','Interpreter','latex','FontSize',ffsize)
box on
set(gca,'TickLabelInterpreter','Latex')


%% gamma OE OPT
% subfigure d
subplot('position', [px(4) py(4) xa ya])
% plot data
semilogy(mu0vec,sqrt(J_vec_xs),stylename{1}), hold on
% highlight values
for jk = 1:length(mu0)
plot(mu0vec(mu0vec == mu0(jk)),sqrt(J_vec_xs(mu0vec == mu0(jk))),stylename2{jk},'markersize',fmsize,'linewidth',flsize);
end

% axis settings
axis([-0.01,0.9,10^0,10^2])
ylabel('$$ \gamma_{OE} (x_{s-opt}) $$','Interpreter','latex','FontSize',ffsize)
text(-0.1826,1e2,'$(d)$','Interpreter','latex','FontSize',ffsize)
xlabel('$$ \mu_0 $$','Interpreter','latex','FontSize',ffsize)
set(gca,'FontSize',ffsize')
set(gca,'TickLabelInterpreter','Latex')
%% gamma FI OPT
% subfigure e
subplot('position', [px(5) py(5) xa ya])
% plot data
semilogy(mu0vec,sqrt(J_vec_xa),stylename{2}), hold on
% highlight values
for jk = 1:length(mu0)
plot(mu0vec(mu0vec == mu0(jk)),sqrt(J_vec_xa(mu0vec == mu0(jk))),stylename2{jk},'markersize',fmsize,'linewidth',flsize);
end

% axis settings
axis([-0.01,0.9,10^0,10^2])
xlabel('$$ \mu_0 $$','Interpreter','latex','FontSize',ffsize), 
ylabel('$$ \gamma_{FI} (x_{a-opt}) $$','Interpreter','latex','FontSize',ffsize)
set(gca,'FontSize',ffsize')
text(-0.1826,1e2,'$(e)$','Interpreter','latex','FontSize',ffsize)
set(gcf,'renderer','painters')
set(gca,'TickLabelInterpreter','Latex')

%% save figure
% print('-depsc2','fig4')