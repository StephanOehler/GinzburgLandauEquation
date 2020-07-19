% This file generates Figure 1 of J. Fluid Mech. (2018), vol. 854, pp.
% 34-55

% clear workspace
clear

% Scale figure based on OS
if ismac
    figure('color','white','position',[0,0,sqrt(2),1/2]*400)
else
    figure('color','white','position',[0,0,sqrt(2),1/2]*600)
end

% margins for figure
xl = 0.075; % left margin
xr = 0.015; % right margin
xg = 0.10; % gap between plots
xa = (1 - xl - xr - sum(xg))/ (length(xg) + 1); % width of each plot
yt = 0.03; % top margin
yb = 0.17; % bottom margin
ya = (1 - yt - yb); % height of each plot

px(1) = xl; % x coordinate of (a)
px(2) = px(1) + xa + xg; % x coordinate of (b)
py([1,2]) = yb; % y coordinate of (a) and (b)


ffsize = 13; % font size of figure
flsize = get(0,'DefaultLineLinewidth')*1; % linewidth of figure
fmsize = get(0,'DefaultLineMarkerSize')*1.33; % marker size of figure



addpath('/Users/oehlers/Google Drive/University/Ph.D.Work/CodeCollection/Workfolder/Matfiles')
gl = CGLe.spectral('SupCrit');

gl.L = 20;
Nmax = 5;

mu0_vec = [0.71,0.56,0.41];

subplot('position', [px(2) py(2) xa ya]), hold on
xmin = -0.8;
xmax = 0.40001;
ymin = -0.8;
ymax = 0;
h1 = plot([xmin,xmax],[0,0],'k:','markersize',fmsize);
h2 = plot([0,0],[ymin,ymax],'k:','markersize',fmsize);
axis tight;
%%
subplot('position', [px(1) py(1) xa ya]), hold on
box on
axis([-20,20,-1,1])
xmin = -20;
xmax = 20;
ymin = -1;
ymax = 1;
h3 = plot([xmin,xmax],[0,0],'k:','markersize',fmsize);
h4 = plot([0,0],[ymin,ymax],'k:','markersize',fmsize);

%%
plotoptions1 = {'k-','k--','k-.'};
plotoptions2 = {'ks','kd','ko'};

for i = 1:length(mu0_vec)
    gl.mu0 = mu0_vec(i);
    %%
    subplot('position', [px(1) py(1) xa ya])
    wimax = gl.mux + gl.cu^2;
    plot(gl.xgrid,wimax,plotoptions1{i})
    if real(gl.x_unstwnmb) ~= 0
        plot(gl.x_unstwnmb * [-1,1],0,plotoptions2{i},'markersize',fmsize)
    end
    %%
    subplot('position', [px(2) py(2) xa ya])
    lambda = gl.mu0 - gl.cu^2 - (gl.nu^2)/4/gl.gam - ( [0:(Nmax - 1)] + 0.5)*gl.hh;
    plot(lambda,plotoptions2{i},'markersize',fmsize)
end
%%

subplot('position', [px(1) py(1) xa ya])
set(gca,'FontSize',ffsize)
xlabel('$x$','interpreter','latex','FontSize',ffsize),
ylabel('$\omega_{i,max}$','interpreter','latex','FontSize',ffsize),
set(gca,'TickLabelInterpreter','Latex')
text(-27.5,1,'$(a)$','interpreter','latex','FontSize',ffsize)
subplot('position', [px(2) py(2) xa ya])
set(gca,'FontSize',ffsize)
xlabel('Re','interpreter','latex','FontSize',ffsize),
ylabel('Im','interpreter','latex','FontSize',ffsize),
text(-1.025,0,'$(b)$','interpreter','latex','FontSize',ffsize)
set(gca,'TickLabelInterpreter','Latex')
box on

print('fig1','-depsc')

