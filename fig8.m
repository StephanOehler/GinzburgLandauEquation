% This file generates Figure 7 of J. Fluid Mech. (2018), vol. 854, pp.
% 34-55

% clear workspace
clear; close all


% margins for figure
xl = 0.055; % left margin
xr = 0.015; % right margin
xg = 0.10; % gap between plots
xa = (1 - xl - xr - sum(xg))/ (length(xg) + 1); % width of each plot
yt = 0.03; % top margin
yb = 0.17; % bottom margin
ya = (1 - yt - yb); % height of each plot

px(1) = xl; % x coordinate of (a)
px(2) = px(1) + xa + xg; % x coordinate of (b)
py([1,2]) = yb;% y coordinate of (a) and (b)


% figure propterties
ffsize = 13; % font size of figure
flsize = get(0,'DefaultLineLinewidth')*1; % linewidth of figure
fmsize = get(0,'DefaultLineMarkerSize')*1.33; % marker size of figure
options = [fmsize,flsize];
markers = {'kx','ko','kv','k^'};
markersacling = [1.5,1,1.1,1.25];
linetype = {'--','-','-.',':'};
markerheight = 0;

gl = CGLe.dynamic();

gl.U = 2; % U: convection constant
gl.cu = 0.2; % cc: most unstable wavenumber
gl.cd = -1; % cd: dispersion
gl.mu0 = 0.41; % mu0: control paramter
gl.mu2 = -0.01; % mu2: non-parallel parameter
gl.nx = 150; % nx: number of Hermite points
% width of actuator and sensor
gl.var_a = 0.4 * sqrt(2);
gl.var_s = 0.4 * sqrt(2);
gl.type = 'non-parallel';

% width of full domain
gl.L = 25;

% weights for batycentric interpolation and domain
b = (bary_weights(gl.nx - 1));
xx = -25:0.001:25;

% upstream and downstream branch
x1=-gl.x_unstwnmb;
x2=gl.x_unstwnmb;


% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,sqrt(2),1/2]*monitor_size(4)/1.5)

%% part a

% create subplot
subplot('position', [px(1) py(1) xa ya]), hold on
set(gca,'FontSize',ffsize)

% vector of sensor locations
xs_vec = [0,2.1053,7.28,gl.x_unstwnmb];

% plot figure a
for i = 1:length(xs_vec)

    % set sensor location
    gl.x_s = xs_vec(i);  
    
    % Solve the continuous algebratic ricatti equation:
    [Y,~,~] = care(gl.A',gl.Cy',gl.Bw*gl.Bw',gl.Vp05^2);
    
    % calculate H2 norm for each gridpoint
    gamma_2 = gl.Q \(diag(gl.Cw * Y * gl.Cw')); % diag(Y);
    
    % do barycentric intertoplation
    ff = bary_interp_new(gl.xgrid , xx , b , sqrt(gamma_2));
    
    
    % plot lines and marker
    plot(xx, ff,['k',linetype{i}])
    plot(xs_vec(i),markerheight,markers{i},'markersize',markersacling(i) * options(1),'linewidth',options(2))
    
end

% apply box, labels, limits and LaTex Style
box on
xlim([-20,20])
xlabel('$x$','interpreter','latex','FontSize',ffsize)
ylabel('$\epsilon_{OE}$','interpreter','latex','FontSize',ffsize)
ylim([0,4.2])
text(-25,4.2,'(a)','interpreter','latex','FontSize',ffsize)
set(gca,'TickLabelInterpreter','Latex')

%% part 2

% vector of sensor locations
xa_vec = [0, -2.02,-7.28,-gl.x_unstwnmb];

% create subplot
subplot('position', [px(2) py(2) xa ya]), hold on
set(gca,'FontSize',ffsize)

for i = 1:length(xa_vec)

    % set actuator location
    gl.x_a = xa_vec(i);

    
    % Solve the continuous algebratic ricatti equation:
    [X,~,~] = care(gl.A,gl.Bu,gl.Cw'*gl.Cw,gl.Rp05^2);
    
    % calculate H2 norm for each gridpoint
    gamma_2 = gl.Q \( diag(gl.Bw' * X * gl.Bw) ); %diag(X);
    
    % do barycentric intertoplation
    ff = bary_interp_new(gl.xgrid , xx , b , sqrt(gamma_2));
    
    % plot lines and marker
    plot(xx,ff,['k',linetype{i}])
    plot(xa_vec(i),markerheight,markers{i},'markersize',markersacling(i) * options(1),'linewidth',options(2))

end
% apply box, labels, limits and LaTex Style
xlim([-20,20])
ylim([0,4.2])
box on
xlabel('$x$','interpreter','latex','FontSize',ffsize)
%ylabel('$\kappa_{FI,d}$','interpreter','latex','FontSize',ffsize)
ylabel('$\epsilon_{FI}$','interpreter','latex','FontSize',ffsize)
text(-25,4.2,'(b)','interpreter','latex','FontSize',ffsize)
set(gca,'TickLabelInterpreter','Latex')
%set(0,'DefaultLineLineWidth',default_LW,'DefaultLineMarkerSize',default_MS,'DefaultAxesFontSize',default_AF);

% save figure
% print('fig8','-depsc')

