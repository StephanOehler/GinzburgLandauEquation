% This file generates the support files for Figure 6 of
% J. Fluid Mech. (2018), vol. 854, pp. 34-55

clear

% add high-level folder to path
addpath('../..')

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

ord = 10;


%%
% delay range
tau_step = 0.25;
tau_min = 0;
tau_max = 6.25;
tau_vec = tau_min:tau_step:tau_max;


%%%% Brute Force Range
%
x_step_dense = 0.5;
x_min = -2.5;
x_max = 2.5;
xa_vecs = x_min:x_step_dense:x_max;
[xa_vecs,~] = meshgrid(xa_vecs,tau_vec);
xs_vecs = x_min:x_step_dense:x_max;
[xs_vecs,~] = meshgrid(xs_vecs,tau_vec);

disp('Run 1')
gamma_2_mat = loop_over_delay(gl,xa_vecs,xs_vecs,tau_vec,ord);
%save('LQG_surf05')

% Update to resolution 0.1
[xa_vecs,xs_vecs] = update_xa_vecs(0.1,0.5,xa_vecs,xs_vecs,tau_vec,gamma_2_mat);


% Run at resolution 0.1
disp('Run 2')
gamma_2_mat = loop_over_delay(gl,xa_vecs,xs_vecs,tau_vec,ord);
%save('LQG_surf01')


% Update to resolution 0.01
[xa_vecs,xs_vecs] = update_xa_vecs(0.01,0.1,xa_vecs,xs_vecs,tau_vec,gamma_2_mat);

% Run at resolution 0.01
disp('Run 3')
gamma_2_mat = loop_over_delay(gl,xa_vecs,xs_vecs,tau_vec,ord);
%save('LQG_surf001')


% Update to resolution 0.001
[xa_vecs,xs_vecs] = update_xa_vecs(0.001,0.01,xa_vecs,xs_vecs,tau_vec,gamma_2_mat);

% Run at resolution 0.001
disp('Run 4')
gamma_2_mat = loop_over_delay(gl,xa_vecs,xs_vecs,tau_vec,ord);
save('LQG_time_delay','gamma_2_mat','tau_vec','xa_vecs','xs_vecs')
%%

% This function updates the vectors that contains the range of steps for a
% particular delay
function [xa_vecs,xs_vecs] = update_xa_vecs(x_step_dense,x_step_prev,xa_vecs,xs_vecs,tau_vec,gamma_2_mat)

% ensure gamma is real and that there are no bad results
gamma_2_mat = real(gamma_2_mat);
gamma_2_mat(gamma_2_mat == 0) = NaN;

% find the minimum H2 norm
[a,b] = min(gamma_2_mat);
a = squeeze(a);
b = squeeze(b);
[~,d] = min(a);
bb = zeros(size(d));
for i = 1:length(tau_vec)
    bb(i) = b(d(i),i);
end

% find the corresponding optimal sensor amd actuator locations
xs_range = zeros(size(d));
for i = 1 : length(d)
    xs_range(i) = xs_vecs(i,d(i));
end
xa_range = zeros(size(bb));
for i = 1 : length(bb)
    xa_range(i) = xa_vecs(i,bb(i));
end


% change the number of steps based on the resulutions
x_points = x_step_prev / x_step_dense * 2 + 1;


% calculate the new sensor and actuator locations that are to be checked for performance
xs_vecs = zeros(length(tau_vec), x_points);
xa_vecs = zeros(length(tau_vec), x_points);
for i = 1:length(tau_vec)
    xs_vecs(i,:) = (xs_range(i) - x_step_prev) : x_step_dense : (xs_range(i) + x_step_prev);
end
for i = 1:length(tau_vec)
    xa_vecs(i,:) = (xa_range(i) - x_step_prev) : x_step_dense : (xa_range(i) + x_step_prev);
end


end






%%

function gamma_2_mat = loop_over_delay(gl,xa_vecs,xs_vecs,tau_vec,ord)

xa_points = size(xa_vecs,2);
xs_points = size(xs_vecs,2);

gamma_2_mat = [xa_points,xs_points,size(xa_vecs,2)];
tic
for i = 1:xa_points
    for j = 1:xs_points
        for k = 1:length(tau_vec)
            
            % update sensor location
            gl.x_a = xa_vecs(k,i);
            gl.x_s = xs_vecs(k,j);
            
            % introduce delay
            sysx = ss(gl.A,[gl.Bu,gl.Bw],[gl.Cy;gl.Cw],0,'OutputDelay',[tau_vec(k),zeros(1,gl.nx)]);
            sysx = pade(sysx,ord);
            
            % find new delayed matrices
            Adel = sysx.A;
            B2del = sysx.B(:,1);
            C2del = sysx.C(1,:);
            B1del = sysx.B(:,2:end);
            C1del = sysx.C(2:end,:);
            
            %%% Solve the continuous algebratic ricatti equation:
            [X,~,~] = care(Adel,B2del,C1del'*C1del,gl.Rp05^2);
            [Y,~,~] = care(Adel',C2del',B1del*B1del',gl.Vp05^2);
            
            %%
            Gamma_2 = trace(C1del * Y * C1del') + trace((gl.Vp05^2) \ C2del * Y*X*Y * C2del');
            gamma_2_mat(i,j,k) = sqrt(real(Gamma_2));
        end
    end
    time_total = (toc / i) * (xa_points - i) ;
    disp(['i = ',num2str(i),'/',num2str(xa_points),', approx. ', num2str(round(time_total)), ' sec, ', num2str(round (time_total/ 60)), ' min, or ', num2str(time_total/ 60 / 60), ' hours left'])
end

end