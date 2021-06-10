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

% Brute force range
x_step = 0.1;
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

for k = 1:length(mu0)
    gl.mu0 = mu0(k); % mu0: update stability paramter
    
    % allocate memory
    gamma_2_IOC = zeros(length(x_vec),length(x_vec));
    
    % simulation timer
    time_total = 0;
    tic
    
    
    for i = 1:x_vec_length
        for j = 1:x_vec_length
            
            % set sensor location
            gl.x_s = x_vec(i);
            gl.x_a = x_vec(j);
            
            % calculate norm
            gamma_2_IOC(i,j) = gl.gammaIO;
            
        end
        % update simulation timer
        if rem(i,5) == 1
            time_total = (toc / i) * (x_vec_length - i) ;
            disp(['Run: ',num2str(k),', i = ',num2str(i),'/',num2str(x_vec_length),', approx. ', num2str(round(time_total)), ' sec, ', num2str(round (time_total/ 60)), ' min, or ', num2str(time_total/ 60 / 60), ' hours left'])
        end
        
    end
    toc
    save(['LQG_',savename{k}],'x_v ec','gamma_2_IOC')
end

warning('on',warningid)