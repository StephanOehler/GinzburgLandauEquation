% This file generates the support files for Figure 6 of
% J. Fluid Mech. (2018), vol. 854, pp. 34-55
clear

% add high-level folder to path
addpath('..')

% Load in default parameters for the GL equation
gl = CGLe.dynamic('SupCrit');

% change limits of the flow
gl.L = 25; % default is 20

% other default parameters are:

% gl.U = 2; % U: convection constant
% gl.cu = 0.2; % cc: most unstable wavenumber
% gl.cd = -1; % cd: dispersion
% gl.mu2 = -0.01; % mu2: non-parallel parameter
% gl.nx = 150; % nx: number of grid points
% gl.var_a = 0.4 * sqrt(2); % width of actuators
% gl.var_s = 0.4 * sqrt(2); % width of sensors
% gl.type = 'non-parallel'; % type of flow
% gl.scaling = 'linear'; % scaling of flow: 


mu0vec = -0.01:0.001:1;
mu0vec_length = length(mu0vec);
x_s_vec = zeros(size(mu0vec));
x_a_vec = zeros(size(mu0vec));
J_vec = zeros(size(mu0vec));

x_s_start = 0; 
x_a_start = 0; 


%  Set options for fminunc (In the paper, we use a diffferent minimisation
%  algorithm. Unfortunately, we cannot share that part of the code and,
%  therefore, quasi-newton is used instead)
options = optimoptions(@fminunc,'Algorithm','quasi-newton', 'MaxIter',500,'Display','off');
% MAKRE SURE Global Optimization Toolbox is installed
tic
for i = 1:length(mu0vec)
    gl.mu0 = mu0vec(i);

%%


%  Run fminunc to obtain the sensor poostion 
%  This function will return theta and the cost 
[x_vec_temp, J_vec(i)] = ...
	fminunc(@(x)(costFunctionIO(x,gl,[0.001,0.001])), [x_s_start,x_a_start], options);

x_s_vec(i) = x_vec_temp(1);
x_a_vec(i) = x_vec_temp(2);

%  speed up search by updating IC
x_s_start = x_s_vec(i);
x_a_start = x_a_vec(i);

%%
   % update simulation timer
    if rem(i,10) == 1
        time_total = (toc / i) * (mu0vec_length - i) ;
        disp(['i = ',num2str(i),'/',num2str(mu0vec_length),', approx. ', num2str(round(time_total)), ' sec, ', num2str(round (time_total/ 60)), ' min, or ', num2str(time_total/ 60 / 60), ' hours left'])
    end
end
%% 
save('LQG_mu0range','mu0vec','x_s_vec','x_a_vec','J_vec')


