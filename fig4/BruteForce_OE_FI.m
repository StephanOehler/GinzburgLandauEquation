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

% Brute force range
x_step = 0.05;
x_min = -15;
x_max = 15;
x_vec = x_min:x_step:x_max;
x_vec_length = length(x_vec);

%%

%supress ricatti stability warning (for extreme cases x_vec ~ |15|)
warningid = 'Control:foundation:RiccatiAccuracy';
warning('off',warningid)

mu0 = [0.41,0.56,0.71];
savename = {'041','056','071'};

for j = 1:length(mu0)
gl.mu0 = mu0(j); % mu0: update stability paramter

% allocate memory
gamma_2_OE = zeros(length(x_vec),1);
gamma_2_FI = gamma_2_OE;

% simulation timer    
time_total = 0;
tic
for i = 1:length(x_vec)    
        
        % set sensor location
        gl.x_s = x_vec(i);
        gl.x_a = x_vec(i);
        
        % calculate norm
        gamma_2_OE(i) = gl.gammaOE;
        gamma_2_FI(i) = gl.gammaFI;
              
        % update simulation timer    
        if rem(i,100) == 1
        time_total = (toc / i) * (x_vec_length - i) ;
        disp(['Run: ',num2str(1),', i = ',num2str(i),'/',num2str(x_vec_length),', approx. ', num2str(round(time_total)), ' sec, ', num2str(round (time_total/ 60)), ' min, or ', num2str(time_total/ 60 / 60), ' hours left'])
        end
end
toc 
save(['LQE_LQR_',savename{j}],'x_vec','gamma_2_OE','gamma_2_FI')
end

warning('on',warningid)