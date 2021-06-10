function [J, grad] = costFunctionIO(x,gl,dx)


x_s = x(1);
x_a = x(2);

dxs = dx(1);
dxa = dx(end);



C2 = exp(-((gl.xgrid-x_s)/gl.var_s).^2).*diag(gl.Q); C2 = C2';
C2dC2 = exp(-((gl.xgrid-( x_s + dxs ))/gl.var_s).^2).*diag(gl.Q); C2dC2 = C2dC2';

% delta
dC2 = C2dC2 - C2;

% actuator vector
B2 = exp(-((gl.xgrid-x_a)/gl.var_a).^2);

% perturbed actuator vector
B2dB2 = exp(-((gl.xgrid-( x_a + dxa ))/gl.var_a).^2);

% delta
dB2 = B2dB2 - B2;


% ricatti equations to get H2 norms
[X,~,~] = care(gl.A,B2,gl.Cw'*gl.Cw,gl.Rp05^2);
[Y,~,~] = care(gl.A',C2',gl.Bw*gl.Bw',gl.Vp05^2);
%J = trace(gl.Cw * Y * gl.Cw') + trace( (gl.Vp05^-2) * C2 * Y*X*Y * C2');
J = real(trace(gl.Bw' * X * gl.Bw) + trace( inv(gl.Rp05^2) * B2' * X*Y*X * B2));
J = real(J);


%% DgamDxa

% calculates the gradient (see paper for more information)
A1 = gl.A - B2 * inv(gl.Rp05^2) * B2' * X;
QQ = dB2 * inv(gl.Rp05^2)* B2' + B2 * inv(gl.Rp05^2) * dB2';
QQ = - X * QQ * X;
delta_x = lyap(A1',QQ);

delta_Gamma_2 = trace( (gl.Vp05^-2) * C2 * Y * delta_x * Y * C2');
grad1 = real(delta_Gamma_2) / dxa;

%% DgamDxs

% calculates the gradient (see paper for more information)
A2 = gl.A - Y * C2' * inv(gl.Vp05^2) * C2;
QQ = dC2' * inv(gl.Vp05^2)* C2 + C2' * inv(gl.Vp05^2) * dC2;
QQ = - Y * QQ * Y;
delta_y = lyap(A2,QQ);

delta_Gamma_2 = trace( inv(gl.Rp05^2) * B2' * X * delta_y * X * B2);
grad2 = real(delta_Gamma_2) / dxs;

%%
grad = [grad1, grad2];
end