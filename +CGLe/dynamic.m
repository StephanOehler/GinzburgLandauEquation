classdef dynamic < CGLe.spectral
    
    % This Class for represents the complex Ginzburg Landau equation
    
    % It is used in J. Fluid Mech. (2018), vol. 854, pp. 34-55
    % See the paper for more information in the parameters and formulas
    
    % Also See J. Fluid Mech. (2011), vol.681, pp. 241-260
    
    properties
        x_a % actuator location
        x_s % sensor location
        var_a % variance of gaussian function representing the actuator
        var_s % variance of gaussian function representing the actuator
        Rp05 = 1/7 % Penalisation actutor force ($\alpha$ in the JFM 2018)
        Vp05 = (1e-3) %covariance sensor noice (Vp05 = $v$ in the JFM 2018)
    end
    
    properties (Dependent, Access = public)
        Bu % actutor matrix
        Cy % sensor matrix
        Bw % actuation everywhere
        Cw % measurements everywhere
        gammaOE % $gamma_{OE}$ is a H2 norm in the JFM 2019
        gammaFI % $gamma_{OE}$ is a H2 norm in the JFM 2019
        gammaIO % $gamma_{OE}$ is a H2 norm in the JFM 2019
    end
    
    properties (Dependent)
        p2 % H2 norm of the uncontrolled system
    end
    
    methods
        %% Setup function - default overwrite
        function theCGLE = dynamic(Setup)
            if nargin == 1
                if strcmp(Setup,'SubCrit')
                    theCGLE.U = 2;
                    theCGLE.cu = 0.2;
                    theCGLE.cd = -1.0;
                    theCGLE.mu0 = 0.38;
                    theCGLE.mu2 = -0.01;
                    theCGLE.nx = 150;
                    theCGLE.L = 20;
                    theCGLE.scaling = 'linear';
                    %%% Actuator and Sensor
                    theCGLE.x_a = -1;
                    theCGLE.var_a = 0.4 * sqrt(2);
                    theCGLE.x_s = -theCGLE.x_a;
                    theCGLE.var_s = theCGLE.var_a;
                end
                if strcmp(Setup,'SupCrit')
                    theCGLE.U = 2;
                    theCGLE.cu = 0.2;
                    theCGLE.cd = -1.0;
                    theCGLE.mu0 = 0.41;
                    theCGLE.mu2 = -0.01;
                    theCGLE.nx = 150;
                    theCGLE.L = 20;
                    theCGLE.scaling = 'linear';
                    %%% Actuator and Sensor
                    theCGLE.x_a = -1;
                    theCGLE.var_a = 0.4 * sqrt(2);
                    theCGLE.x_s = -theCGLE.x_a;
                    theCGLE.var_s = theCGLE.var_a;
                end
                if strcmp(Setup,'Stable')
                    theCGLE.U = 2;
                    theCGLE.cu = 0.2;
                    theCGLE.cd = -1.0;
                    theCGLE.mu0 = 0.41;
                    theCGLE.mu2 = -0.01;
                    theCGLE.nx = 150;
                    theCGLE.L = 20;
                    theCGLE.scaling = 'linear';
                    %%% Actuator and Sensor
                    theCGLE.x_a = -1;
                    theCGLE.var_a = 0.4 * sqrt(2);
                    theCGLE.x_s = -theCGLE.x_a;
                    theCGLE.var_s = theCGLE.var_a;
                end
            end
        end
        
        function identify(theCGLE)
            disp(['I am a dynamical complex Ginzburg-Landau system discretised with nx = ',num2str(theCGLE.nx),' chebshev points'])
        end
        
        %% Set Bu (actuator matrix)
        function Bu = get.Bu(theCGLE)
            if theCGLE.var_a ~= 0
                Bu = zeros(length(theCGLE.xgrid), length(theCGLE.x_a));
                nx = length(theCGLE.x_a);
                nv = length(theCGLE.var_a);
                % loop through each actuator location
                for i = 1:nx
                    j = rem(i,nv); if j == 0, j = nv; end
                    Bu(:,i) = exp( -((theCGLE.xgrid - theCGLE.x_a(i)) ...
                        / theCGLE.var_a(j)) .^ 2);
                end
            elseif theCGLE.var_a == 0
                Bu = eye(size(theCGLE.A));
            end
        end
        %% Set C (sensor matrix)
        function Cy = get.Cy(theCGLE)
            if theCGLE.var_s ~= 0
                Cy = zeros(length(theCGLE.xgrid), length(theCGLE.x_s));
                nx = length(theCGLE.x_s);
                nv = length(theCGLE.var_s);
                % loop through each sensor location
                for i = 1:nx
                    Cy(:,i) = exp( -((theCGLE.xgrid - theCGLE.x_s(i)) ...
                        / theCGLE.var_s(rem(nx,nv) + 1)) .^ 2);
                end
                Cy = theCGLE.Q * exp( -((theCGLE.xgrid - theCGLE.x_s) ...
                    / theCGLE.var_s ) .^2 );
                Cy = Cy';
            elseif theCGLE.var_s == 0
                Cy = eye(size(theCGLE.A));
                Cy = Cy';
            end
        end
        %% set Cw (Measuring everywhere)
        function Cw = get.Cw(theCGLE)
            Cw = sqrt(theCGLE.Q);
        end
        %% set Bw (Actuation everywhere
        function Bw = get.Bw(theCGLE)
            Bw = sqrt(inv(theCGLE.Q));
        end
        %% Reshape x_a
        function theCGLE = set.x_a(theCGLE,x_a)
            theCGLE.x_a = reshape(x_a,[1,numel(x_a)]);
        end
        %% Reshape x_s
        function theCGLE = set.x_s(theCGLE,x_s)
            theCGLE.x_s = reshape(x_s,[1,numel(x_s)]);
        end
        %% Reshape var_a
        function theCGLE = set.var_a(theCGLE,var_a)
            if  all(var_a > 0) && ~isempty(var_a)
                theCGLE.var_a = reshape(var_a,[1,numel(var_a)]);
            elseif length(var_a) == 1 && (var_a == 0)
                theCGLE.var_a = 0;
            else
                error('Positive actuator variance required. When var_a = 0, then an identical force is applied at each grid point in the flow.')
            end
        end
        %% Reshape var_s
        function theCGLE = set.var_s(theCGLE,var_s)
            if  all(var_s > 0) && ~isempty(var_s)
                theCGLE.var_s = reshape(var_s,[1,numel(var_s)]);
            elseif length(var_s) == 1 && (var_s == 0)
                theCGLE.var_s = 0;
            else
                error('Positive actuator variance required. When var_a = 0, then an identical force is applied at each grid point in the flow.')
            end
        end
        %% Restricted Rp05
        function theCGLE = set.Rp05(theCGLE,Rp05)
            if  any(Rp05 > 0)
                theCGLE.Rp05 = Rp05;
            else
                error('R cannot be negative')
            end
        end
        %% Restricted Vp05
        function theCGLE = set.Vp05(theCGLE,Vp05)
            if  any(Vp05 > 0)
                theCGLE.Vp05 = Vp05;
            else
                error('Vp05 cannot be negative')
            end
        end
        %%
        function impulse(ThisSystem,tvm,nt)
            % dt - time step
            % tvm - max time
            % tvm - time vector, if dt not given
            if nargin == 2
                impulse(ThisSystem.P,tvm)
            elseif nargin == 3
                impulse(ThisSystem.P,linspace(0,tvm,nt))
            else
                impulse(ThisSystem.P)
            end
        end
        
        function bode(ThisSystem)
            bode(ThisSystem.P)
        end
        function pzplot(ThisSystem)
            pzplot(ThisSystem.P)
        end
        %% Get uncontrolled model
        function p2 = get.p2(theCGLE)
            Y = lyap(theCGLE.A,theCGLE.Bw * theCGLE.Bw');
            p2 = sqrt(trace(theCGLE.Cw * Y * theCGLE.Cw'));
        end
        %% Set gammaIO
        function gammaIO = get.gammaIO(theCGLE)
            B1 = [theCGLE.Bw, zeros(theCGLE.nx,1)];
            C1 = [theCGLE.Cw; zeros(1,theCGLE.nx)];
            
            [X,~,~] = care(theCGLE.A,theCGLE.Bu,C1'*C1,theCGLE.Rp05^2);
            [Y,~,~] = care(theCGLE.A',theCGLE.Cy',B1*B1',theCGLE.Vp05^2);
            
            gammaIO = sqrt(trace(C1 * Y * C1') + trace((theCGLE.Vp05^2)...
                \ theCGLE.Cy * Y*X*Y * theCGLE.Cy'));
        end
        %% Set gammaOE
        function gammaOE = get.gammaOE(theCGLE)
            B1 = [theCGLE.Bw, zeros(theCGLE.nx,1)];
            C1 =  theCGLE.Cw;
            [Y,~,~] = care(theCGLE.A',theCGLE.Cy',B1*B1',theCGLE.Vp05^2);
            gammaOE = sqrt(trace(C1 * Y * C1'));
        end
        %% Set gammaFI
        function gammaFI = get.gammaFI(theCGLE)
            C1 = theCGLE.Cw;
            B1 = theCGLE.Bw;
            [X,~,~] = care(theCGLE.A,theCGLE.Bu,C1'*C1,theCGLE.Rp05^2);
            gammaFI = sqrt(trace(B1' * X * B1));
        end
        
        %% Get most unstable mode
        function S_min = Get_RHP_p_z_S_min(theCGLE)
            [z, p] = zpkdata(theCGLE.P,'v');
            
            z_vec = z(real(z) > 0);
            p_vec = p(real(p) > 0);
            if isempty(z_vec) || isempty(p_vec)
                S_min = NaN;
            else
                Cj = ones(1,length(z_vec));
                for j = 1:length(z_vec)
                    for k = 1:length(p_vec)
                        Cj(j) = Cj(j) * abs( z_vec(j) + conj(p_vec(k))) ...
                            / abs(z_vec(j) - p_vec(k) );
                    end
                end
                S_min = max(Cj);
            end
        end
    end
end
