classdef basic
    % This Class for represents the complex Ginzburg Landau equation
    
    % It is used in J. Fluid Mech. (2018), vol. 854, pp. 34-55
    % See the paper for more information in the parameters
    
    % The complex Ginzburg-Landau equation needs to know its:
    % U: convection constant
    % cu: most unstable wavenumber
    % cd: dispersion
    % mu0: control paramter
    % mu2: non-parallel parameter
    % nx: number of spatial points
    
    %%  Private Properties
    properties
        U
        cu
        cd
        mu0
        mu2
    end
    %% Dependent properties
    properties (Dependent, Access = public)
        %%  Extra terms
        nu                  % nu - first order scaling term
        gam                 % gamma - second order scaling term
        Umax                % Overall velocity
        chi                 % decay rate of global scalings
        hh                  % h constant, part of eigenvalue expression
        %% Criteria for absolute instability and global instability
        mut                 % Locally absolute unstability criterion
        muc                 % Globally unstable criterion
        %% Branches
        % Up and downstream branches based on App. Mech. Rev. (2009), Vol
        % 62 / 020803-1, equation(7) and following paragraph
        x_unstregn
        % Up and downstream branches based on App. Mech. Rev. (2009), Vol
        % 62 / 020803-1, equation(3) and following paragraph
        % (branches used in J. Fluid Mech. (2018), vol. 854, pp. 34-55)
        x_unstwnmb
        %%
    end
    %% Methods
    methods
        %% Setup function - default overwrite
        function thisCGLE = basic(Setup)
            if nargin == 1
                % Values are taken from App. Mech. Rev. (2009), Vol 62 / 
                % 020803-1
                if strcmp(Setup,'SupCrit')
                    thisCGLE.U = 2;
                    thisCGLE.cu = 0.2;
                    thisCGLE.cd = -1.0;
                    thisCGLE.mu0 = 0.41;
                    thisCGLE.mu2 = -0.01;
                end
            end
        end
        %% Identify
        function identify(thisCGLE)
            disp(['I am a complex Ginzburg-Landau system with U_{max} = ',num2str(thisCGLE.Umax),' L/T'])
        end
        %% Get gamma
        function gam = get.gam(theCGLE)
            gam = 1 + 1i * theCGLE.cd;
        end
        %% Get chi
        function chi = get.chi(theCGLE)
            chi = (-theCGLE.mu2 / (2 * theCGLE.gam))^(0.25);
        end
        %% Get nu
        function nu = get.nu(theCGLE)
            nu = theCGLE.U + 2 * theCGLE.cu * real(theCGLE.gam) * 1j;
        end
        %% Get hh
        function hh = get.hh(theCGLE)
            hh = sqrt(- 2 * theCGLE.mu2 * theCGLE.gam);
        end
        %% Get Umax
        function Umax = get.Umax(theCGLE)
            Umax = theCGLE.U + 2 * theCGLE.cu * theCGLE.cd;
        end
        %% Get mut
        function mut = get.mut(theCGLE)
            mut = (theCGLE.Umax)^2 / (4*abs(theCGLE.gam)^2);
        end
        %% Get muc
        function muc = get.muc(theCGLE)
            muc = theCGLE.mut + abs(theCGLE.hh)/2*cos(angle(theCGLE.gam)/2);
        end
        %% Get x_unstregn
        function x_unstregn = get.x_unstregn(theCGLE)
               x_unstregn = sqrt( - 2 * (theCGLE.mu0 - theCGLE.cu ^ 2 ) / theCGLE.mu2);
        end
        %% Get x_unstwnmb
        function x_unstwnmb = get.x_unstwnmb(theCGLE)
            x_unstwnmb = sqrt( - 2 * theCGLE.mu0 / theCGLE.mu2);
        end        
        %% Restrict mu2
        function theCGLE = set.mu2(theCGLE,mu2)
            if mu2 < 0
                theCGLE.mu2 = mu2;
            else
                error('mu2 must be < 0')
            end
        end
    end
end