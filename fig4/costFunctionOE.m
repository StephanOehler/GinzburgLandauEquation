function [J, grad] = costFunctionOE(xs,gl,dxs)
% calculate gradient  (J. Fluid Mech. (2018), vol. 854, pp. 34-55)

% sensor vector
C2 = exp(-((gl.xgrid-xs)/gl.var_s).^2).*diag(gl.Q); C2 = C2';

% perturbed sensor vector
C2dC2 =  exp(-((gl.xgrid-( xs + dxs ))/gl.var_s).^2).*diag(gl.Q); C2dC2 = C2dC2';

% delta
dC2 = C2dC2 - C2;

% ricatti equation to get H2 norms
[Y,~,~] = care(gl.A',C2',gl.Bw*gl.Bw',gl.Vp05^2);
J = trace(gl.Cw * Y * (gl.Cw'));
J = real(J);


%% DgamDxs

% calculates the gradient (see paper for more information)
A2 = gl.A - Y * C2' * inv(gl.Vp05^2) * C2;
QQ = dC2' * inv(gl.Vp05^2)* C2 + C2' * inv(gl.Vp05^2) * dC2;
QQ = - Y * QQ * Y;
delta_y = lyap(A2,QQ);

delta_Gamma_2 = trace( gl.Cw * delta_y * gl.Cw' );
grad = real(delta_Gamma_2) / dxs;
end