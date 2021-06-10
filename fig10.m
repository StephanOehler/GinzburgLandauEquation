% This file generates Figure 1 of J. Fluid Mech. (2018), vol. 854, pp.
% 34-55

% clear workspace
clear; close all

%% setup figure

% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,sqrt(2),1/2]*monitor_size(4)/1.5)


% margins for figure
xl = 0.07; % left margin
xr = 0.015; % right margin
xg = 0.09; % gap between plots
xa = (1 - xl - xr - sum(xg))/ (length(xg) + 1); % width of each plot
yt = 0.03; % top margin
yb = 0.18; % bottom margin
ya = (1 - yt - yb); % height of each plot


% x coordinates for subfigures
px(1) = xl; % x coordinate of (a)
px(2) = px(1) + xa + xg; % x coordinate of (b)
py([1,2]) = yb; % y coordinate of (a) and (b)

%%

Umax = 1.6;
%%
cut = 26;

ffsize = 13;
%%
load('fig10/LQR_delay')

[a,b] = min(gamma_2_mat,[],2);
tau_vec_length = length(tau_vec);
xa_opt = zeros(1,tau_vec_length);
for i = 1:tau_vec_length
 xa_opt(i) = xa_vecs(i,b(i));
end

%%
subplot('Position',[px(1),py(1),xa,ya]); hold on
plot(xa_opt,tau_vec * Umax ,'k--')
%%

subplot('Position',[px(2),py(2),xa,ya]);
semilogy(tau_vec * Umax  ,a,'k--'), hold on
%%



%% 
load('fig10/LQE_delay')

[a,b] = min(gamma_2_mat,[],2);
tau_vec_length = length(tau_vec);
xs_opt = size(tau_vec);
for i = 1:tau_vec_length
 xs_opt(i) = xs_vecs(i,b(i));
end
%%
subplot('Position',[px(1),py(1),xa,ya]);
plot(xs_opt,tau_vec * Umax ,'k-')
%%
subplot('Position',[px(2),py(2),xa,ya]);
semilogy(tau_vec * Umax, (a(1:cut)),'k-')

%%

load('fig10/LQG_time_delay')


xs_vec = zeros(1,tau_vec_length);
for i = 1:tau_vec_length
 xs_vec(i) = xs_vecs(i,b(i));
end
xa_vec = zeros(1,tau_vec_length);
for i = 1:tau_vec_length
 xa_vec(i) = xa_vecs(i,b(i));
end

%%
subplot('Position',[px(1),py(1),xa,ya]);
plot(xs_vec(1:cut),tau_vec(1:cut) * 1.6 ,'k-.')
plot(xa_vec(1:cut),tau_vec(1:cut) * 1.6 ,'k:')
set(gca,'TickLabelInterpreter','Latex')

%%
subplot('Position',[px(2),py(2),xa,ya]);
c = min(gamma_2_mat,[],2);
d = min(gamma_2_mat,[],1);
semilogy(tau_vec * 1.6  ,squeeze(d),'k-.')



%% labels

subplot('Position',[px(1),py(1),xa,ya]);
set(gca,'FontSize',ffsize)
box(gca,'on')
text(-4.0,10,'$(a)$','interpreter','latex','FontSize',ffsize)
xlabel('$x$','interpreter','latex'),
ylabel(' $\tau U_{max}$','interpreter','latex')
xlim([-3,3])

subplot('Position',[px(2),py(2),xa,ya]);
set(gca,'FontSize',ffsize)
ylabel('$\gamma$','interpreter','latex'), xlabel(' $\tau U_{max}$','interpreter','latex')
set(gca,'YTick',[5,10,15,20,25])
ylim([5,25])

text(-5/3,25,'$(b)$','interpreter','latex','FontSize',ffsize)

set(gca,'TickLabelInterpreter','Latex')

 print('fig10','-depsc2')

