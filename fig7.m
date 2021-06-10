% This file generates Figure 7 of J. Fluid Mech. (2018), vol. 854, pp.
% 34-55

% clear workspace
clear; close all


% margins for figure
xl = 0.01; % left margin
xr = 0.015; % right margin
xa = (1 - xl - xr); % width of each plot
yt = 0.03; % top margin
yb = 0.18; % bottom margin
ya = (1 - yt - yb); % height of each plot


ffsize = 13; % font size of figure
flsize = 1*get(0,'DefaultLineLinewidth'); % linewidth of figure
fmsize = 2*get(0,'DefaultLineMarkerSize'); % marker size of figure
%%
addpath('external_files')
gl = CGLe.basic('SupCrit');

%%

% range of domain
x = -18:0.01:18;

% number of grid points
nx = length(x);

%% calculate eigenmode and adjoint eigenmode
phi = abs(exp( gl.nu / 2 / gl.gam * x - gl.chi^2 * x.^2 / 2 ));
phi_adj = abs(conj(phi) .* exp( - conj(gl.nu)  .* x ./ conj(gl.gam) ));

% normalise
phi = phi / max(phi);
phi_adj = phi_adj / max(phi_adj);

% calculate wavemaker and normalise
wavemaker = phi .* phi_adj;
wavemaker = wavemaker / max(wavemaker);


% stability branch points
x1 = -gl.x_unstwnmb;
x2 = gl.x_unstwnmb;


%% Plot global modes

% open figure and set plot options
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,sqrt(2),1/2]*monitor_size(4)/1.5)
subplot('Position',[xl,yb,xa,ya]),hold on
set(gca,'FontSize',ffsize)
linestyle = {'k-','k--','k-.','k:'};
options = [fmsize,flsize];
markerheight = -0.1;


% calculate grey background
y1 = 0;
y2 = 1;
fill([x1 x2 x2 x1 x1], [y1 y1 y2 y2 y1],[1 1 1] * 0.95,'edgecolor',[0.9 0.9 0.9])
plot([x1 x2 x2 x1 x1], [y1 y1 y2 y2 y1],'k','LineWidth',0.1,'LineStyle',':')

%n plot lines
plot(x,phi,linestyle{1})
plot(x,phi_adj,linestyle{2})
plot(x,wavemaker,linestyle{3},'linewidth',0.5)

% find peaks
[a,b] = max(phi); [c,d] = max(phi_adj); 

% plot peaks
plot(x([b,d]),markerheight,'kv','markersize',1.1 * options(1),'linewidth',options(2))

% plot branch points
plot([x1,x2],markerheight,'k^','markersize',1.25*options(1),'linewidth',options(2))
plot(0,markerheight,'kx','markersize',1.5*options(1),'linewidth',options(2))
plot([-1.1,1.1],markerheight,'ko','markersize',options(1),'linewidth',options(2))

% plot branch labels
text(x1 - 0.6,-0.3,'$$ X_{I} $$','Interpreter','Latex','FontSize',ffsize)
text(x2 - 0.6,-0.3,'$$ X_{II} $$','Interpreter','Latex','FontSize',ffsize)

% plot x-axis label
xlabel('Spatial Position','FontSize',ffsize,'Interpreter','Latex')
set(gca, 'box','off','YTickLabel',[],'YTick',[],'YColor', 'white')

% set limits
ylim([-0.1,1])
xlim([min(x), max(x)])

% switch to LaTex style
set(gca,'TickLabelInterpreter','Latex')

% print figure
% print('fig7','-depsc')