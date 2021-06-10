function [J, grad] = costFunctionFI(x_a,gl,dxa)
% calculate gradient  (J. Fluid Mech. (2018), vol. 854, pp. 34-55)

% actuator vector
B2 = exp(-((gl.xgrid-x_a)/gl.var_a).^2);

% perturbed actuator vector
B2dB2 = exp(-((gl.xgrid-( x_a + dxa ))/gl.var_a).^2);

% delta
dB2 = B2dB2 - B2;

% ricatti equation to get H2 norms
[X,~,~] = care(gl.A,B2,gl.Cw'*gl.Cw,gl.Rp05^2); 
J = trace(gl.Bw' * X * gl.Bw);
J = real(J);


%% DgamDxa

% calculates the gradient (see paper for more information)
A1 = gl.A - B2 * inv(gl.Rp05^2) * B2' * X;
QQ = dB2 * inv(gl.Rp05^2)* B2' + B2 * inv(gl.Rp05^2) * dB2';
QQ = - X * QQ * X;
delta_x = lyap(A1',QQ);

delta_Gamma_2 = trace(gl.Bw' * delta_x * gl.Bw);
grad = real(delta_Gamma_2) / dxa;
end