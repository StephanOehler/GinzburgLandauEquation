classdef spectral < CGLe.basic
    % This Class for represents the complex Ginzburg Landau equation
    
    % It is used in J. Fluid Mech. (2018), vol. 854, pp. 34-55
    % See the paper for more information in the parameters and formulas
    
    % The complex Ginzburg-Landau equation needs to know its:
    
    % U: convection constant
    % cu: most unstable wavenumber
    % cd: dispersion
    % mu0: control paramter
    % mu2: non-parallel parameter
    % nx: number of spatial points in the chebyshev grid with boundary
    % conditions applied
    % L: linear scaling factor of the chebyshev grid
    
    % scaling: scaling type: chose {'linear','fraction','atan','seminf'}
    % (see book:
    %  Stability and Transition in Shear Flows
    %  by PJ Schmid and DS Henningson in Appendix A and
    %  Chebyshev and Fourier Spectral Methods
    %  JP Boyd. U Michigan, 2000. Chapter 16.4 p 326-327)
    
    % Type: 'non-parallel' and 'parallel'
    % (see App. Mech. Rev. (2009), Vol 62 / 020803-1)
    
    
    %%  Private Properties
    properties
        %% Discretisation
        nx
        L
        scaling = 'linear';
        type = 'non-parallel';
    end
    %% Dependent properties
    properties (Dependent, Access = public)
        %% Discretised System terms
        xh % unscaled chebyshev grid
        D1 % First order differentiation matrix
        D2 % Second order differentiation matrix
        xgrid % scaled chebyshev grid
        Q % Integration matrix
        %% Result
        mux % quadratic function (see equation (7) in App. Mech. Rev.
        % (2009), Vol 62 / 020803-1)
        omega_imax % The growth rate of the most unstable wavenumber
        A % System matrix A
    end
    %% Unscaled Chebyshev terms
    properties (Dependent, Access = protected)
        xC % unscaled full chebyshev grid
        DC % unscaled full differentiation matrix
        wC % unscaled full clenshaw curtis weights
    end
    %% Methods
    methods
        %% Setup function - default overwrite
        function thisCGLE = spectral(Setup)
            if nargin == 1
                if strcmp(Setup,'SupCrit')
                    
                    % Values are taken from App. Mech. Rev. (2009), Vol 62 /
                    % 020803-1
                    
                    thisCGLE.U = 2;
                    thisCGLE.cu = 0.2;
                    thisCGLE.cd = -1.0;
                    thisCGLE.mu0 = 0.41;
                    thisCGLE.mu2 = -0.01;
                    thisCGLE.nx = 200;
                    thisCGLE.L = 20;
                    thisCGLE.scaling = 'linear';
                end
            end
        end
        %% Identify
        function identify(thisCGLE)
            disp(['I am a complex Ginzburg-Landau system discretised with nx = ',num2str(thisCGLE.nx),' chebyshev points'])
        end
        %% Get xC
        function xC = get.xC(theCGLE)
            xC = gl_chebx(theCGLE.nx + 1);
        end
        %% Get DC
        function DC = get.DC(theCGLE)
            DC = gl_chebD(theCGLE.nx + 1,theCGLE.xC);
        end
        %% Get wC
        function wC = get.wC(theCGLE)
            wC = gl_chebw(theCGLE.nx + 1);
        end
        %% Get xh
        function xh = get.xh(theCGLE)
            xh = theCGLE.xC(2:theCGLE.nx + 1);
        end
        %% Get D1
        function D1 = get.D1(theCGLE)
            D1 = gl_scaleD1(theCGLE.xh,theCGLE.DC,theCGLE.nx,theCGLE.L,theCGLE.scaling);
        end
        %% Get D2
        function D2 = get.D2(theCGLE)
            D2 = gl_scaleD2(theCGLE.xh,theCGLE.DC,theCGLE.nx,theCGLE.L,theCGLE.scaling);
        end
        %% Get xgrid
        function xgrid = get.xgrid(theCGLE)
            xgrid = gl_scale_xgrid(theCGLE.xh,theCGLE.L,theCGLE.scaling);
        end
        %% Get Q (diag(w))
        function Q = get.Q(theCGLE)
            Q = diag(gl_scale_w(theCGLE.xh,theCGLE.wC,theCGLE.L,theCGLE.scaling));
            %Q = diag( ([ diff(theCGLE.xgrid); 0 ] + [ 0; diff(theCGLE.xgrid) ]) / 2 ); % Trapezoidal rule
        end
        %% Get mux
        function mux = get.mux(theCGLE)
            if strcmp(theCGLE.type,{'non-parallel'})
                mux = theCGLE.mu0 - real( theCGLE.gam ) * theCGLE.cu ^ 2 + theCGLE.mu2 * ( theCGLE.xgrid ) .^ 2 / 2 ;
            elseif strcmp(theCGLE.type,{'parallel'})
                mux = theCGLE.mu0 - real( theCGLE.gam ) * theCGLE.cu ^ 2;
            end
        end
        %% Get omegaimax
        function omega_imax = get.omega_imax(theCGLE)
            if strcmp(theCGLE.type,{'non-parallel'})
                omega_imax = theCGLE.mu0 + theCGLE.mu2 * ( theCGLE.xgrid ) .^ 2 / 2 ;
            elseif strcmp(theCGLE.type,{'parallel'})
                omega_imax = theCGLE.mu0;
            end
        end
        %% Get A
        function A = get.A(theCGLE)
            A = - theCGLE.nu * theCGLE.D1 + theCGLE.gam * theCGLE.D2 + diag( theCGLE.mux );
        end
        %% Restrict scaling
        function theCGLE = set.scaling(theCGLE,scaling)
            if  any(strcmp(scaling,{'linear','fraction','atan','semiinf'}))
                theCGLE.scaling = scaling;
            else
                error('scaling must be either: ''linear'',''fraction'',''atan'',''semiinf'' ')
            end
        end
        %% Restrict equation type
        function theCGLE = set.type(theCGLE,type)
            if  any(strcmp(type,{'non-parallel','parallel'}))
                theCGLE.type = type;
            else
                error('type must be either: ''non-parallel'',''parallel'' ')
            end
            if strcmp(type,'parallel')
                warning('This gives the parallel comlex-Ginzburg-Lanadu equation. Results might not be what you expect them to be.')
            end
        end
        %% Restrict nx
        function theCGLE = set.nx(theCGLE,nx)
            if nx > 0
                theCGLE.nx = nx;
            else
                error('nx must be > 0')
            end
        end
    end
end

%%

% See Trefethen, Lloyd N. Spectral methods in MATLAB.
% Society for industrial and applied mathematics, 2000.
% Page 54, code cheb.m

function x = gl_chebx(N)
x = cos( pi * (0:N) / N)';
end

function D = gl_chebD(N,x)
c = [2; ones(N - 1 , 1); 2] .* (-1).^(0:N)';
X = repmat(x , 1 , N + 1 );
dX = X - X';
D = ( c * (1./ c)') ./ (dX + (eye (N + 1)));
D = D - diag( sum(D'));
end


% See Trefethen, Lloyd N. Spectral methods in MATLAB.
% Society for industrial and applied mathematics, 2000.
% Page 54, code clencurt.m


function w = gl_chebw(N)
theta = pi*(0:N)'/N;
w = zeros(1,N+1);
ii = 2:N;
v = ones(N-1,1);
if mod(N,2)==0
    w(1) = 1/(N^2-1);
    w(N+1) = w(1);
    for k=1:N/2-1,
        v = v - 2*cos(2*k*theta( ii))/(4*k^2-1);
    end
    v = v - cos(N*theta(ii))/(N^2-1);
else
    w(1) = 1/N^2; w(N+1) = w(1);
    for k=1:(N-1)/2,
        v = v - 2*cos(2*k*theta(ii))/(4*k^2-1);
    end
end
w(ii) = 2*v/N;
end


%%

% See:
%  Stability and Transition in Shear Flows
%  by PJ Schmid and DS Henningson in Appendix A and
%  Chebyshev and Fourier Spectral Methods
%  JP Boyd. U Michigan, 2000. Chapter 16.4 p 326-327)
%  and the pdf: GL_spectral_code

function D1 = gl_scaleD1(xh,DM,nx,L,mode)
D1 = DM(2:nx + 1, 2:nx + 1); % First order differential
if strcmp(mode,'linear')
    D1 = D1 / L;
elseif strcmp(mode,'fraction')
    D1 = 1 / L *  D1 * diag( (1 - xh.^2).^ (3/2) );
elseif strcmp(mode,'atan')
    D1 = D1 * diag( 1 -  xh.^2 ) / L;
elseif strcmp(mode,'semiinf')
    D1 = D1 / 2 / L * diag((1 - xh).^2);
else
    D1 = [];
    
end
end

function D2 = gl_scaleD2(xh,DM,nx,L,mode)
D1 = DM(2:nx + 1, 2:nx + 1); % First order differential
D2 = DM^2; % First order differential\
D2 = D2(2:nx + 1, 2:nx + 1);
if strcmp(mode,'linear')
    D2 = D2 / L^2;
elseif strcmp(mode,'fraction')
    D2 = 1 / L^2 * ( -D1 * diag (-3 * xh .* sqrt( 1 - xh.^2 ).^2) + D2 * diag( (1 - xh.^2).^3 ) );
elseif strcmp(mode,'atan')
    D2 = 1  / L^2 * ( -D1 * diag (  2 * xh .* (1 - xh.^2)  ) + D2 * diag( (1 - xh.^2).^2 ) );
elseif strcmp(mode,'semiinf')
    D2 = 1 / L^2 * ( -D1 * diag ( xh ) + D2 * diag( (1 - xh.^2) ));
    
else
    D2 = [];
end
end


function xgrid = gl_scale_xgrid(xh,L,mode)
if strcmp(mode,'linear')
    xgrid = L * xh;
elseif strcmp(mode,'fraction')
    xgrid =  L * xh ./ sqrt(1 - xh.^2);
elseif strcmp(mode,'atan')
    xgrid = L * atanh(xh);
elseif strcmp(mode,'semiinf')
    xgrid = L * (1 + xh) ./ (1 - xh);
else
    xgrid = [];
end
end


function w = gl_scale_w(xh,w,L,mode)
if strcmp(mode,'linear')
    w = L * w(2:end-1)';
elseif strcmp(mode,'fraction')
    w = L ./ ((1 - xh.^2).^(3/2)) .* w(2:end-1)';
elseif strcmp(mode,'atan')
    w = L ./ (1 - xh.^2) .*  w(2:end-1)';
elseif strcmp(mode,'semiinf')
    w = L *  2  ./ (1 - xh).^2 .* w(2:end-1)';
    
else
    w = [];
end

%%
end