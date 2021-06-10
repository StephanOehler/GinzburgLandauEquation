% This file generates Figure 1 of J. Fluid Mech. (2018), vol. 854, pp.
% 34-55

% clear workspace
clear; close all

% set width of domain
xlmin = -20;
xumax = 20;

% set resolution of domain
xres = 0.00001;

% get each point in the domain
x = xlmin:xres:xumax;

% Load in default parameters
gl = CGLe.basic('SupCrit');

% Extract default parameters
cd = gl.cd;
U = gl.U;
mu0 = gl.mu0;
mu2 = gl.mu2;
cc = gl.cd;
gam =gl.gam;


% calculate the first four eigenmodes (from the analytical hermite
% polynomial based solution)
phi = [];
phi(:,1) = abs(exp( gl.nu / 2 / gl.gam * x - gl.chi^2 * x.^2 / 2 ));
phi(:,2) = abs(exp( gl.nu / 2 / gl.gam * x - gl.chi^2 * x.^2 / 2 ) .* 2 .* gl.chi .* x );
phi(:,3) = abs(exp( gl.nu / 2 / gl.gam * x - gl.chi^2 * x.^2 / 2 ) .* (4 .* gl.chi .^ 2 .* x .^ 2 - 2) );
phi(:,4) = abs(exp( gl.nu / 2 / gl.gam * x - gl.chi^2 * x.^2 / 2 )) .* (8 .* gl.chi .^ 3 .* x .^ 3 - 12 .* gl.chi .* x);

% calculate the first four adjoint eigenmodes (from the analytical hermite
% polynomial based solution)
phi_adj = [];
phi_adj(:,1) = conj(phi(:,1))' .* exp( - conj(gl.nu)  .* x ./ conj(gl.gam) );
phi_adj(:,2) = conj(phi(:,2))' .* exp( - conj(gl.nu)  .* x ./ conj(gl.gam) );
phi_adj(:,3) = conj(phi(:,3))' .* exp( - conj(gl.nu)  .* x ./ conj(gl.gam) );
phi_adj(:,4) = conj(phi(:,4))' .* exp( - conj(gl.nu)  .* x ./ conj(gl.gam) );
phi_adj = abs(phi_adj);

%% setup figure

% Scale figure based on OS

% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,sqrt(2),3/4]*monitor_size(4)/1.5)

% margins for figure
xl = 0.085; % left margin
xr = 0.015; % right margin
xg = 0.10; % gap between plots
xa = (1 - xl - xr - sum(xg))/ (length(xg) + 1); % width of each plot
yt = 0.02; % top margin
yb = 0.115; % bottom margin
yg = 0.05;
ya = (1 - yt - yb - sum(yg)) / (length(yg) + 1); % height of each plot



px(1) = xl; % x coordinate of (a)
px(2) = px(1) + xa + xg; % x coordinate of (b)
px(3:4) = px(1:2); % x coordinates of (c,d)

py(3) = yb; % y coordinate of (c)
py(1) = py(3) + ya + yg; % y coordinate of (a)
py([2;4]) = py([1,3]); % y coordinate of (b,d)


ffsize = 13; % font size of figure
flsize = get(0,'DefaultLineLinewidth')*1; % linewidth of figure


%% Plot figure

linestyle = {'k-','k--','k-.',':'};

%% Plot eigenmodes

for j = [1,3]
    
    % call figure
    subplot('position', [px(j) py(j) xa ya])
    
    % set the style of the lines
    
    % loop through and plot the first MM eigenmodes
    for MM = 1:3
        hold on;
        plot(x,(abs(phi(:,MM))/max(abs(phi(:,MM)))),linestyle{MM})
    end
    
    % Scale axes, label axes, label figure:

    set(gca,'box','on')
    set(gca,'FontSize',ffsize)
    ylabel('$$\phi$$','Interpreter','Latex','FontSize',ffsize)
    if j == 1
        axis([xlmin+5,xumax,0,1])
        set(gca,'yscale','linear')
        text(-22,1,'$(a)$','Interpreter','Latex','FontSize',ffsize)
        set(gca,'XTick',(-20:5:20))
        set(gca,'XTickLabel',{'','','','','','','','',''})
    else
        set(gca,'yscale','log')
        axis([xlmin+5,xumax,10^-6,1])
        text(-22,1,'$(c)$','Interpreter','Latex','FontSize',ffsize)
        set(gca,'XTick',(-20:5:20))
        xlabel('$x$','Interpreter','Latex','FontSize',ffsize)
        set(gca,'YTick',[1e-6,1e-3,1])
    end
    set(gca,'TickLabelInterpreter','Latex')
end


%% Plot adjoint eigenmodes


for j = [2,4]
    
    % call figure
    subplot('position', [px(j) py(j) xa ya])
    
    % loop through and plot the first MM adjoint eigenmodes
    for MM = 1:3
        hold on;
        plot(x,(abs(phi_adj(:,MM))/max(abs(phi_adj(:,MM)))),linestyle{MM})%,xh,real(Us(:,MM))/max(abs(Us(:,MM))),'k--',xh,-abs(Us(:,MM))/max(abs(Us(:,MM))),'k','LineWidth',1);axis([xh(end) -xh(end) -1.2 1.2 ]);
    end
    
    % Scale axes, label axes, label figure:
    set(gca,'box','on')
    set(gca,'FontSize',ffsize)
    ylabel('$$\psi$$','Interpreter','Latex','FontSize',ffsize)
    if j == 2
        axis([xlmin,xumax-5,0,1])
        set(gca,'yscale','linear')
        text(-27,1,'$(b)$','Interpreter','Latex','FontSize',ffsize)
        set(gca,'XTick',(-20:5:20))
        set(gca,'XTickLabel',{'','','','','','','','',''})
    else
        set(gca,'yscale','log')
        axis([xlmin,xumax-5,10^-6,1])
        text(-27,1,'$(d)$','Interpreter','Latex','FontSize',ffsize)
        set(gca,'XTick',(-20:5:20))
        xlabel('$x$','Interpreter','Latex','FontSize',ffsize)
        set(gca,'YTick',[1e-6,1e-3,1])
    end
    set(gca,'TickLabelInterpreter','Latex')
    
end

%% print esp file
% print('fig2','-depsc')