% This file generates the support files for Figure 9 of
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

%%

mu0_vec = 0:0.001:0.9;

Gamma_IO = NaN(size(mu0_vec));
Gamma_FI = Gamma_IO;
Gamma_OE = Gamma_IO;

gl.x_a = 0;%-1.12;
gl.x_s = 0;%1.09;

for i = 1:length(mu0_vec)
    
    
    gl.mu0 = mu0_vec(i);
    
    Gamma_IO(i) = gl.gammaIO;
    Gamma_FI(i) = gl.gammaFI;
    Gamma_OE(i) = gl.gammaOE;
    
end
save('Simple_h2_norm_range')