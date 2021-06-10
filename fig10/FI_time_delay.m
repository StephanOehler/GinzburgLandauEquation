% This file generates the support files for Figure 6 of
% J. Fluid Mech. (2018), vol. 854, pp. 34-55

clear

% add high-level folder to path
addpath('../..')

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


% order of pade time delay
ord = 10;


%%
% delay range
tau_step = 0.25;
tau_min = 0;
tau_max = 6.25;
tau_vec = tau_min:tau_step:tau_max;

% Brute Force Range
x_step_dense = 0.25;
x_min = -3;
x_max = 3;
xa_vecs = x_min:x_step_dense:x_max;
[xa_vecs,~] = meshgrid(xa_vecs,tau_vec);


% Run at resolution 0.25
disp('Run 1')
gamma_2_mat = loop_over_delay(gl,xa_vecs,tau_vec,ord);
%save('LQE_surf025')

% Update to resolution 0.1
xa_vecs = update_xa_vecs(0.1,x_step_dense,xa_vecs,gamma_2_mat);

% Run at resolution 0.1
disp('Run 2')
gamma_2_mat = loop_over_delay(gl,xa_vecs,tau_vec,ord);
%save('LQE_surf01')
%%

% Update to resolution 0.01
xa_vecs = update_xa_vecs(0.01,0.1,xa_vecs,gamma_2_mat);

% Run at resolution 0.01
disp('Run 3')
gamma_2_mat = loop_over_delay(gl,xa_vecs,tau_vec,ord);
%save('LQE_surf001')
%%

% Update to resolution 0.001
xa_vecs = update_xa_vecs(0.001,0.01,xa_vecs,gamma_2_mat);

% Run at resolution 0.001
disp('Run 4')
gamma_2_mat = loop_over_delay(gl,xa_vecs,tau_vec,ord);
save('LQR_delay','gamma_2_mat','tau_vec','xa_vecs')

%%
% This function updates the vectors that contains the range of steps for a
% particular delay
function [xa_vecs] = update_xa_vecs(x_step_dense,x_step_prev,xa_vecs,gamma_2_mat)

% ensure gamma is real and that there are no bad results
gamma_2_mat = real(gamma_2_mat);
gamma_2_mat(gamma_2_mat == 0) = NaN;

% find the minimum H2 norm
[~,d] = min(log10(gamma_2_mat'));

% find the corresponding optimal actuator location
xa_range = zeros(size(d));
for i = 1 : length(d)
    xa_range(i) = xa_vecs(i,d(i));
end

% change the number of steps based on the resulutions
x_points = x_step_prev / x_step_dense * 2 + 1;

% calculate the new actuator locations that are to be checked for performance
xa_vecs = zeros(size(xa_vecs,1), x_points);
for i = 1:size(xa_vecs,1)
    xa_vecs(i,:) = (xa_range(i) - x_step_prev) : x_step_dense : (xa_range(i) + x_step_prev);
end
end


function gamma_2_mat = loop_over_delay(gl,xa_vecs,tau_vec,ord)

gamma_2_mat = size(xa_vecs);
tic
for i = 1:size(xa_vecs,2)
    for j = 1:length(tau_vec)
        
        % update sensor location
        gl.x_a = xa_vecs(j,i);
        
        % introduce delay
        sysx = ss(gl.A,[gl.Bu, gl.Bw],gl.Cw,0,'InputDelay',[tau_vec(j),zeros(1,gl.nx)]);
        sysx = pade(sysx,ord);
        
        % find new delayed matrices
        Atdl = sysx.A;
        B2tdl = sysx.B(:,1);
        B1tdl = sysx.B(:,2:end);
        C1tdl = sysx.C;
        
        % calculate the new gamma value
        [X,~,~] = care(Atdl,B2tdl,C1tdl'*C1tdl,gl.Rp05^2);
        
        % calculate the H2 norm
        Gamma_2 = trace(B1tdl' * X * B1tdl);
        
        % ssave the H2 norm
        gamma_2_mat(j,i) = sqrt(Gamma_2);
        
    end
    if rem(i,5) == 1
        time_total = (toc / i) * (size(xa_vecs,2) - i) ;
        disp(['i = ',num2str(i),'/',num2str(size(xa_vecs,2)),', approx. ', num2str(round(time_total)), ' sec, ', num2str(round (time_total/ 60)), ' min, or ', num2str(time_total/ 60 / 60), ' hours left'])
    end
end

end