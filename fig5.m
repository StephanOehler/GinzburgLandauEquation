% This file generates Figure 5 of J. Fluid Mech. (2018), vol. 854, pp.
% 34-55

% clear workspace
clear; close all

% add support files
addpath('external_files')

% create figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,sqrt(2),1/2]*monitor_size(4)/1.5)

% margins for figure
xl = 0.08; % left margin
xr = 0.065; % right margin
xg = [0.04 0.01]; % gap between plots
xc = 0.03; % coloarbar
xa = (1 - xl - xr - xc - sum(xg))/(length(xg)+0); % width of each plot
yt = 0.03; % top margin
yb = 0.18; % bottom margin
ya = (1 - yt - yb); % height of each plot
px(1) = xl;
px(2) = px(1) + xa + xg(1);
px(3) = px(2) + xa + xg(2);
py([1,2,3]) = yb;


ffsize = 13; % font size of figure
flsize = get(0,'DefaultLineLinewidth')*1; % linewidth of figure
fmsize = get(0,'DefaultLineMarkerSize')*1.33; % marker size of figure

%
load fig5/LQELQR_RMS.mat
X_II_vec = real(X_II_vec);
%% OE

% find the peak (once for the left and once for the right half)
[idxx,~] = size(OE_2_mat_bary);
[~,b] = max(OE_2_mat_bary(1:floor(idxx/2),:));
[c,d] = max(OE_2_mat_bary(ceil(idxx/2+1000):idxx,:));
d = d + floor(idxx/2+1000);


% plot lines and surface
h1 = subplot('position', [px(1) py(1) xa ya]);
p1 = pcolor( xx(1001:end-1000),mu0vec,log10(sqrt(real(OE_2_mat_bary(1001:end-1000,:)')) ));  shading flat, hold on
plot(xx(b),mu0vec,'k--')
plot(smooth(c,xx(d),.5,'rloess'),mu0vec,'k--')
caxis([-1,1.5])
set(gca,'Fontsize',ffsize)
set(h1,'Box','On','linewidth',1, 'TickDir', 'Out','fontsize',ffsize)
set(p1,'edgecolor','none')


% plot branch
plot(X_II_vec,mu0vec,'w:','linewidth',1.25)
plot(-X_II_vec,mu0vec,'w:','linewidth',1.25)

% axis settings
xlabel('$$ x $$','interpreter','latex','Fontsize',ffsize)
ylabel('$$ \mu_0 $$','interpreter','latex','Fontsize',ffsize)
text(-27.5,0.9,'$(a)$','interpreter','latex','Fontsize',ffsize)
ylim([-0.01,0.9])
set(gca,'TickLabelInterpreter','Latex')

%% FI

% find the peak (once for the left and once for the right half)
[a,b] = max(FI_2_mat_bary(1:floor(idxx/2-1000),:));
[c,d] = max(FI_2_mat_bary(ceil(idxx/2):idxx,:));
d = d + floor(idxx/2);

% plot lines and surface
h2 = subplot('position', [px(2) py(2) xa ya]);
p2 = pcolor( xx(1001:end-1000),mu0vec,log10(sqrt(real(FI_2_mat_bary(1001:end-1000,:)')) )); shading flat, hold on
plot(smooth(a,xx(b),.5,'rloess'),mu0vec,'k--')
plot(xx(d),mu0vec,'k--')
set(h2,'Box','On','linewidth',1, 'TickDir', 'Out','fontsize',ffsize)
set(p2,'edgecolor','none')

caxis([-1,1.5])
set(gca,'Fontsize',ffsize)


% COLORBAR
colormap(brewermap([],'YlGnBu'))
cbh = colorbar;
set(cbh, 'Position', [px(3) py(1) xc ya])
set(cbh,'Ticks',[-1,-0.5,0,0.5,1,1.5],'Ticklength',0.02)

% branches
plot(X_II_vec,mu0vec,'w:','linewidth',1.25)
plot(-X_II_vec,mu0vec,'w:','linewidth',1.25)

% axis settings
xlabel('$$ x $$','interpreter','latex','Fontsize',ffsize)
set(gca,'YTickLabel',[])
for ij = 1:length(cbh.TickLabels)
    cbh.TickLabels{ij} = strcat('10^{',cbh.TickLabels{ij},'}');
end
ylim([-0.01,0.9])
text(-23.5,0.9,'$(b)$','interpreter','latex','Fontsize',ffsize)
set(gca,'TickLabelInterpreter','Latex')

%%
% print('fig5','-depsc')
