% This file generates the support files for Figure 4 of
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



mu0vec = -0.01:0.005:1;
mu0vec_length = length(mu0vec);
x_s_vec = zeros(size(mu0vec));
J_vec_xs = zeros(size(mu0vec));

x_a_vec = zeros(size(mu0vec));
J_vec_xa = zeros(size(mu0vec));


x_s_start = 0; 
x_a_start = 0; 


%  Set options for fminunc (In the paper, we use a diffferent minimisation
%  algorithm. Unfortunately, we cannot share that part of the code and,
%  therefore, quasi-newton is used instead)
options = optimoptions(@fminunc,'Algorithm','quasi-newton', 'MaxIter', 10000,'Display','none');

tic
for i = 1:length(mu0vec)
    gl.mu0 = mu0vec(i);

%%
% preallocate variables for fast calculation for OE
Parameters.A = gl.A;
Parameters.xgrid = gl.xgrid;
Parameters.Bw = gl.Bw;
Parameters.Cw = gl.Cw;
Parameters.Vp05 = gl.Vp05;
Parameters.var_s = gl.var_s;
Parameters.Q = gl.Q;


%  Run fminunc to obtain the sensor poostion 
%  This function will return theta and the cost 
[x_s_vec(i), J_vec_xs(i)] = ...
	fminunc(@(x_s)(costFunctionOE(x_s,gl,0.001)), x_s_start, options);

%  speed up search by updating IC
x_s_start = x_s_vec(i);


%%
% additional parameters for FI
Parameters.var_a = gl.var_a;
Parameters.Rp05 = gl.Rp05;


%  Run fminunc to obtain the sensor poostion 
%  This function will return theta and the cost 
[x_a_vec(i), J_vec_xa(i)] = ...
	fminunc(@(x_a)(costFunctionFI(x_a,gl,0.001)), x_a_start, options);

%  speed up search by updating IC
x_a_start = x_a_vec(i);


%%
   % update simulation timer
    if rem(i,10) == 1
        time_total = (toc / i) * (mu0vec_length - i) ;
        disp(['i = ',num2str(i),'/',num2str(mu0vec_length),', approx. ', num2str(round(time_total)), ' sec, ', num2str(round (time_total/ 60)), ' min, or ', num2str(time_total/ 60 / 60), ' hours left'])
    end
end
%%
save('LQE_LQR_mu0range','mu0vec','x_s_vec','J_vec_xs','x_a_vec','J_vec_xa')


