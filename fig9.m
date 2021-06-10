% This file generates Figure 6 of J. Fluid Mech. (2018), vol. 854, pp.
% 34-55

% clear workspace
clear; close all

% add support files
addpath('external_files')

% margins for figure
xl = 0.08; % left margin
xr = 0.02; % right margin
xa = (1 - xl - xr); % width of each plot
yt = 0.035; % top margin
yb = 0.185; % bottom margin
ya = (1 - yt - yb); % height of each plot

% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,sqrt(2),1/2]*monitor_size(4)/2) 
subplot('Position',[xl,yb,xa,ya])


ffsize = 13; % font size of figure
flsize = 1*get(0,'DefaultLineLinewidth'); % linewidth of figure
fmsize = 2*get(0,'DefaultLineMarkerSize'); % marker size of figure


% plot stability branches
semilogy([0.395,0.395,NaN,0.555,0.555,NaN,0.705,0.705,NaN,0.865,0.865], ...
  [10^0,10^10,NaN,10^0,10^10,NaN,10^0,10^10,NaN,10^0,10^10],'k:'), hold on;

% load placement at centre fo flow
load fig9/Simple_h2_norm_range.mat

% plot data
gamma_2 = real(Gamma_IO);
hold on
semilogy(mu0_vec,gamma_2,'k-')
set(gca,'FontSize',ffsize)

% label figure
xlabel('$$ \mu_0 $$','Interpreter','Latex')
axis([0.3,0.9,1,1e8])
set(gca,'YTick',[1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8])
ylabel('$$\gamma_{IO}$$','interpreter','latex')

%ylh = get(gca,'ylabel');
%gyl = get(ylh);  
%ylp = get(ylh, 'Position');
%set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
%save2pdf('Discussion1b',gcf,300)

% load optimal lcation data
load fig6/LQG_mu0range.mat

gamma_2 = sqrt(J_vec);

% plot data
hold on,
semilogy(mu0vec,gamma_2,'k--')
set(gca,'TickLabelInterpreter','Latex')

% save figure
print('-depsc2','fig9')
