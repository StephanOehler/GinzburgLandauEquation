%%

addpath('..')
load ../fig4/LQE_LQR_mu0range.mat


% Load in default parameters for the GL equation
gl = CGLe.dynamic('SupCrit');

% change limits of the flow
gl.L = 25; % default is 20

% set up interpolation
OE_2_mat = zeros(gl.nx,length(mu0vec));
a = bary_weights(gl.nx + 1);
xx = gl.L * linspace(-1,+1,10001);
OE_2_mat_bary = zeros(length(xx),length(mu0vec));
X_II_vec = zeros(size(mu0vec));


disp('Generating part a')
for i = 1:length(mu0vec)
    gl.mu0 = mu0vec(i);
    gl.x_s = x_s_vec(i);
    X_II_vec(i) = gl.x_unstwnmb;
    [Y,~,~] = care(gl.A',gl.Cy',gl.Bw*gl.Bw',gl.Vp05^2);
    OE_2_mat(:,i) = gl.Q \ diag( gl.Cw * Y * gl.Cw' );
    OE_2_mat_bary(:,i) = bary_interp(gl.xgrid,xx,a,OE_2_mat(:,i));
end


% set up interpolation
FI_2_mat = zeros(gl.nx,length(mu0vec));
a = bary_weights(gl.nx + 1);
xx = gl.L * linspace(-1,+1,10001);
FI_2_mat_bary = zeros(length(xx),length(mu0vec));
X_II_vec = zeros(size(mu0vec));


disp('Generating part b')
for i = 1:length(mu0vec)
    gl.mu0 = mu0vec(i);
    gl.x_a = x_a_vec(i);
    X_II_vec(i) = gl.x_unstwnmb;
    [X,~,~] = care(gl.A,gl.Bu,gl.Cw'*gl.Cw,gl.Rp05^2); 
    FI_2_mat(:,i) = gl.Q \ diag(gl.Bw' * X * gl.Bw);
    FI_2_mat_bary(:,i) = bary_interp(gl.xgrid,xx,a,FI_2_mat(:,i));
end


%% save data
save('LQELQR_RMS')