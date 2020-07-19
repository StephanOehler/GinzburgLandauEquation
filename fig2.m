clear all; close all; clc

%%%%% REWRITE THIS CODE TO GET THE EIGENMODES VIA HERMITE ANALYTICAL
addpath('/Users/oehlers/Google Drive/University/Ph.D.Work/CodeCollection/Workfolder/Matfiles')

gl = CGLe.basic('SupCrit');


x = -20:0.00001:20;
%x = -20:0.00001:20;


xlmin = min(x);
xumax = max(x);

%x = -15:0.01:15;
nx = length(x);

cd = gl.cd;
U = gl.U;
mu0 = gl.mu0;
mu2 = gl.mu2;
cc = gl.cd;
gam =gl.gam;

phi(:,1) = abs(exp( gl.nu / 2 / gl.gam * x - gl.chi^2 * x.^2 / 2 )); %.* (32*xgrid.^5 - 160 * xgrid.^3 + 120)) ;
phi(:,2) = abs(exp( gl.nu / 2 / gl.gam * x - gl.chi^2 * x.^2 / 2 ) .* 2 .* gl.chi .* x );
phi(:,3) = abs(exp( gl.nu / 2 / gl.gam * x - gl.chi^2 * x.^2 / 2 ) .* (4 .* gl.chi .^ 2 .* x .^ 2 - 2) );
phi(:,4) = abs(exp( gl.nu / 2 / gl.gam * x - gl.chi^2 * x.^2 / 2 )) .* (8 .* gl.chi .^ 3 .* x .^ 3 - 12 .* gl.chi .* x);
%phi_adj = flipud(phi);

phi_adj = [];
phi_adj(:,1) = conj(phi(:,1))' .* exp( - conj(gl.nu)  .* x ./ conj(gl.gam) );
phi_adj(:,2) = conj(phi(:,2))' .* exp( - conj(gl.nu)  .* x ./ conj(gl.gam) );
phi_adj(:,3) = conj(phi(:,3))' .* exp( - conj(gl.nu)  .* x ./ conj(gl.gam) );
phi_adj(:,4) = conj(phi(:,4))' .* exp( - conj(gl.nu)  .* x ./ conj(gl.gam) );
phi_adj = abs(phi_adj);


%
% Unstable domain: branch I and II
%
x1=-sqrt(-2*(mu0-cc^2)/mu2);
x2=sqrt(-2*(mu0-cc^2)/mu2);


%%

%
% Plot global modes
%
figure('color','white','position',[0,0,sqrt(2),3/4]*400), 


    % FIGURE 1 Surface vel mes
% margins
xl = 0.085; % left margin
xr = 0.015; % right margin
xg = 0.10; % gap between plots
xa = (1 - xl - xr - sum(xg))/ (length(xg) + 1); % width of each plot
yt = 0.02; % top margin
yb = 0.115; % bottom margin
yg = 0.05;
ya = (1 - yt - yb - sum(yg)) / (length(yg) + 1); % height of each plot



px(1) = xl;
px(2) = px(1) + xa + xg;
px(3:4) = px(1:2);

py(3) = yb;
py(1) = py(3) + ya + yg;
py([2;4]) = py([1,3]);


ffsize = 13;
flsize = get(0,'DefaultLineLinewidth')*1;


%%




for j = [1,3]
subplot('position', [px(j) py(j) xa ya])


linestyle = {'k-','k--','k-.',':'};
for MM = 1:3
hold on;
plot(x,(abs(phi(:,MM))/max(abs(phi(:,MM)))),linestyle{MM})%,xh,real(Us(:,MM))/max(abs(Us(:,MM))),'k--',xh,-abs(Us(:,MM))/max(abs(Us(:,MM))),'k','LineWidth',1);axis([xh(end) -xh(end) -1.2 1.2 ]);
end
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

%set(gca,'YTick',[])
end
%%
linestyle = {'k-','k--','k-.','-'};

for j = [2,4]
subplot('position', [px(j) py(j) xa ya])
for MM = 1:3
hold on;
plot(x,(abs(phi_adj(:,MM))/max(abs(phi_adj(:,MM)))),linestyle{MM})%,xh,real(Us(:,MM))/max(abs(Us(:,MM))),'k--',xh,-abs(Us(:,MM))/max(abs(Us(:,MM))),'k','LineWidth',1);axis([xh(end) -xh(end) -1.2 1.2 ]);
end
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
%save2pdf('Result4',gcf,300)
print('fig2','-depsc')