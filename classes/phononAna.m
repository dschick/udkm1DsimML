%% phononAna
% The phononAna class simulates phonon dynamics on a 1D sample structure.
%
% Copyright (c) 2013, Daniel Schick, André Bojahr, Marc Herzog, Roman Shayduk, Clemens von Korff Schmising
% All rights reserved.
%
% License: BSD (use/copy/change/redistribute on own risk, mention the authors)

%% Classdef
% Each phononAna instance and all inherited class objects are inherted from 
% the phonon class which follows handle semantics. 
% Hence a copy of such object will not copy the object itself,
% but only a handle to that object.
classdef phononAna < phonon
    %% Properties
    properties (SetAccess=public,GetAccess=public)
    end%properties
    %% Methods
    methods
        %% Constructor
        % Is executed each time an instance of this class is created. Only
        % the _structure_ input is obligatory.
        function obj = phononAna(structure,varargin)
            obj = obj@phonon(structure,varargin{:});
        end%function
        
        %% Display
        % This method is called to display informations of the instance.
        function disp(obj)
            disp('Analytical phonon simulation properties:');
            % call the parent display method
            disp@phonon(obj);
        end%function        
                
        %% getStrainMap
        % Returns a _strainMap_ for the sample structure. If no _strainMap_
        % is saved it is caluclated.
        function [strainMap X V A B sticksSubSystems] = getStrainMap(obj,time,tempMap,deltaTempMap)
            % create a unique hash
            hash = obj.getHash(time,tempMap,deltaTempMap);
            % create the file name to look for
            filename = fullfile(obj.cacheDir, ['strainMapAna_' hash '.mat']);
            if exist(filename,'file') && ~obj.forceRecalc
                % file exists so load it 
                load(filename);
                obj.dispMessage(['_strainMap_ loaded from file ' filename]);
            else
                % file does not exist so calculate and save
                [strainMap X V A B sticksSubSystems] = obj.calcStrainMap(time,tempMap,deltaTempMap);
                save(filename, 'strainMap', 'X', 'V', 'A', 'B', 'sticksSubSystems');
                obj.dispMessage(['_strainMap_ saved to file ' filename]);
            end%if
        end%function
        
        %% calcStrainMap
        % Calculates the _strainMap_ of the sample structure for a given 
        % _tempMap_ and _deltaTempMap_ and _time_ vector. Further details 
        % are given in Ref. [1]. Within the linear chain of $N$ masses 
        % ($m_i$) at position $z_i$ coupled with spring constants $k_i$ one 
        % can formulate the differential equation of motion as follow:
        %
        % $$ m_i\ddot{x}_i = -k_i(x_i-x_{i-1})-k_{i+1}(x_i-x_{i+1}) +
        % F_i^{heat}(t) $$
        %
        % Since we only consider nearest-neighbor interaction one can
        % write:
        %
        % $$ \ddot{x}_i = \sum_{n=1}^N \kappa_{i,n} x_n = \Delta_i(t) $$
        %
        % Here $x_i(t) = z_i(t)-z_i^0$ is the shift of each unit cell,
        % $F_i^{heat}(t)$ is the external force (thermal stress) of each 
        % unit cell and $\kappa_{i,i} = -(k_i + k_{i+1})/m_i$, and 
        % $\kappa_{i,i+1} = \kappa_{i+1,i} = k_{i+1}/m_i$.
        %
        % $k_i = m_i\, v_i^2/c_i^2$ is the spring constant and $c_i$ and
        % $v_i$ are the lattice $c$-axis and longitudinal sound velocity of
        % each unit cell respectively.
        % One can rewrite the homogeneous differential equation in matrix 
        % form to obtain the general solution
        %
        % $$ \frac{d^2}{dt^2} X = K\, X $$
        %
        % Here $X = (x_1 \ldots x_N)$ and $K$ is the
        % tri-diagonal matrix of $\kappa$ which is real and symmetric.
        % The differential equation can be solved with the ansatz:
        %
        % $$ X(t) = \sum_j \Xi_j \, (A_j \cos(\omega_j \, t) + B_j \sin(\omega_j \, t)) $$
        %
        % where $\Xi_j = (\xi_1^j \ldots \xi_N^j)$ are the eigenvectors of
        % the matrix $K$. Thus by solving the Eigenproblem for $K$ one 
        % gets the eigenvecotrs $\Xi_j$ and the eigenfrequencies $\omega_j$.
        % From the initial conditions
        %
        % $$ X(0) = \sum_j \Xi_j \, A_j = \Xi \, A \qquad V(0) = \dot{X}(0) 
        %         = \sum_j \Xi_j \, \omega_j\, B_j = \Xi \, \omega \, B $$
        %
        % one can determine the real coefficient vecots $A$ and $B$ in 
        % order to calculate $X(t)$ and $V(t)$ using the ansatz:
        %
        % $$ A = \Xi \setminus X(0) \qquad B = (\Xi \setminus V(0)) ./ \omega $$
        %
        % The external force is implemented as spacer sticks which are
        % inserted into the springs and hence the unit cells have a new 
        % equillibrium  positions $z_i(\infty) = z_i^\infty$. Thus we can 
        % do a coordination transformation:
        %
        % $$ z_i(t) = z_i^0 + x_i(t) = z_i^\infty + x_i^\infty(t) $$
        %
        % and
        %
        % $$ x_i^\infty(t) = z_i^0 - z_i^\infty + x_i(t) $$
        %
        % with the initial condition $x_i(0) = 0$ the becomes
        %
        % $$ x_i^\infty(0) = z_i^0 - z_i^\infty = \sum_{j = i+1}^N l_j $$
        %
        % $x_i^\infty(0)$ is the new initial condition after the excitation
        % where $l_i$ is the length of the i-th spacer stick. The spacer 
        % sticks are calculated from the temperature change and the linear 
        % thermal expansion coefficients. 
        % The actual strain $\epsilon_i(t)$ of each unit cell is calculates 
        % as follows:
        % 
        % $$ \epsilon_i(t) = [ \Delta x_i(t) + l_i) ] / c_i $$
        %
        % with $\Delta x_i = x_i - x_{i-1}$. The stick $l_i$ have to be
        % added here, because $x_i$ has been transformed into the new
        % coordinate system $x_i^\infty$.
        function [strainMap X V A B sticksSubSystems] = calcStrainMap(obj,time,tempMap,deltaTempMap)
            tic
            % initialize
            N = obj.S.getNumberOfUnitCells; % nb of unit cells
            M = length(time);               % nb of time steps
            
            time0       = time(1); % initial time
            cAxises     = obj.S.getUnitCellPropertyVector('cAxis');
            X           = zeros(M,N); % shifts of the unitCells
            V           = zeros(M,N); % velocities of the unitCells
            A           = zeros(M,N); % coefficient vector for eigenwert solution
            B           = zeros(M,N); % coefficient vector for eigenwert solution
            strainMap   = zeros(M,N); % the restulting strain pattern of the unitCells
            
            % check tempMaps
            [tempMap, deltaTempMap] = obj.checkTempMaps(tempMap,deltaTempMap,time);
            
            % calculate the sticks due to heat expansion first for all time
            % steps
            obj.dispMessage('Calculating linear thermal expansion ...');
            [sticks, sticksSubSystems] = obj.calcSticksFromTempMap(tempMap,deltaTempMap);
                        
            if obj.onlyheat
                % no coherent dynamics so calculate the strain directly
                strainMap = sticks./repmat(cAxises',size(sticks,1),1);
            else
                % solve the eigenproblem for the structure to obtains the
                % eigenvectors Xi and eigenfreqeuencies omega for the N 
                % coupled differential equations
                [Xi, omega] = obj.solveEigenproblem();
            
                % calculate the actual strain pattern with the solution of the
                % eigenproblem and the external force (sticks, thermal stress) 
                obj.dispMessage('Calculating _strainMap_ ...');
                obj.progressBar('Please wait... ');
                % traverse time
                for i=1:M
                    obj.progressBar(i/M*100); % plot the progress

                    dt = time(i)-time0; % this is the time step

                    % calculate the current shift X and velocity V of all 
                    % unitCells using the ansatz
                    X(i,:)  = Xi*(         A(i,:)'.*cos(omega*dt) + B(i,:)'.*sin(omega*dt));
                    V(i,:)  = Xi*(omega.*(-A(i,:)'.*sin(omega*dt) + B(i,:)'.*cos(omega*dt)));
                    % remember the velocities and shifts as ic for the next
                    % time step
                    X0      = X(i,:)';
                    V0      = V(i,:)';

                    % the strain can only be calculated for N-1 unitCells, so
                    % we neglect the last one
                    if i > 1
                        strainMap(i,1:N-1) = (diff(X(i,:),1,2)+sticks(i-1,1:N-1))./cAxises(1:N-1)';
                    else
                        % initial sticks are zero
                        strainMap(i,1:N-1) =  diff(X(i,:),1,2)./cAxises(1:N-1)';
                    end%if

                    % calculate everything for the next step
                    if i < M % check, if there is a next step
                        if find(deltaTempMap(i,:)) % there is a temperature change
                            time0 = time(i); % set new initial time                       
                            
                            % determining the shifts due to inserted sticks
                            % as new ininital conditions
                            if i > 1
                                temp = flipud(cumsum(flipud(sticks(i,:)'-sticks(i-1,:)')));
                            else
                                % initial sticks are zero
                                temp = flipud(cumsum(flipud(sticks(i,:)')));
                            end%if
                            X0 = X0 + vertcat(temp(2:end),0);

                            % determining the cofficient vectors A and B of
                            % the general solution of X(t) using the inital
                            % conditions X0 and V0
                            A(i+1,:) = ( Xi\X0);
                            B(i+1,:) = ((Xi\V0)./omega)';
                        else
                            % no temperature change, so keep the current As,
                            % Bs, and sticks
                            A(i+1,:) = A(i,:);
                            B(i+1,:) = B(i,:);
                        end%if
                    end%if
                end%for
                obj.progressBar('');
            end
            obj.dispMessage('Elapsed time for _strainMap_:',toc);
        end%function
        
        %% solveEigenproblem
        % Creates the real and symmetric $K$ matrix ($N \times N$) of 
        % spring constants $k_i$ and masses $m_i$ and calculates the 
        % eigenvectors $\Xi_j$ and eigenfrequencies $\omega_j$ for the 
        % matrix which are used to calculate the _strainMap_ of the
        % structure.
        % If the result has been save to file, load it from there.
        function [Xi,omega] = solveEigenproblem(obj)
            % create the file name to look for
            filename = fullfile(obj.cacheDir, ['eigenValues_' obj.S.getHash('phonon') '.mat']);
            if exist(filename,'file') && ~obj.forceRecalc
                % file exists so load it 
                load(filename);
                obj.dispMessage(['_eigenValues_ loaded from file ' filename]);
            else
                % no file - so lets calculate everything
                tic
                obj.dispMessage('Calculating _eigenValues_ ...');
                % initialize
                N       = obj.S.getNumberOfUnitCells; % nb of unit cells
                K       = zeros(N,N); %Initializing three-diagonal springs-masses matrix.
                omega   = zeros(N,1); %Initializing a vector for eigenfrequencies
                
                masses       = obj.S.getUnitCellPropertyVector('mass'); % get masses vector
                springConsts = obj.S.getUnitCellPropertyVector('springConst'); % get the first order springs vector
                springConsts = vertcat(0, springConsts(:,1)); % set the first spring free
               
                for i=1:N %Defining main diagonal.
                    K(i,i)=-(springConsts(i)+springConsts(i+1))/masses(i);
                end%for

                for i=2:N %Defining the two other diagonals. Nearest neightbour interaction.
                    K(i,i-1) = springConsts(i)/masses(i);
                    K(i-1,i) = springConsts(i)/masses(i-1);
                end%for
                
                % Determining the eigenvectors and the eigenvalues
                [Xi,lambda] = eig(K);
    
                for i=1:N % calculate the eigenfrequencies from the eigenvalues
                    omega(i)=sqrt(-lambda(i,i));
                end%for
                
                obj.dispMessage('Elapsed time for _eigenValues_:',toc);
                % save the result to file
                save(filename,'Xi', 'omega');
                obj.dispMessage(['_eigenValues_ saved to file ' filename]);
            end%if
        end%function
                
        %% getEnergyPerEigenmode
        % Returns the energy per Eigenmode of the coherent phonons of
        % the 1D sample sorted and unsorted.
        %
        % $$ E_j = \frac{1}{2} (A^2_j + B^2_j)\, \omega_j^2\, m_j \, \| \Xi_j\|^2 $$
        %
        % Frequencies are in [Hz] and energy per mode in [J].
        function [omegaSort ESort omega E] = getEnergyPerEigenmode(obj,A,B)
            % initialize
            N       = obj.S.getNumberOfUnitCells; % nb of unit cells
            M       = size(A,1); % nb of time steps
            E       = zeros(M,N);
            ESort   = zeros(M,N);
            masses  = obj.S.getUnitCellPropertyVector('mass'); % mass vector of unitCells
            
            % get the eigenVectors and eigenFrequencies
            [Xi,omega] = obj.solveEigenproblem(); 
            
            % sort the frequencies and remeber the permutation of indicies
            [omegaSort sortIndex] = sort(omega);            
            
            % traverse time
            for i=1:M
                % calculate the energy for the jth mode
                E(i,:) = 0.5 * (A(i,:)'.^2 + B(i,:)'.^2).* omega(:).^2.*masses(:) .* sum(Xi.^2,1)';
                % sort the energies according to the frequencies
                ESort(i,:) = E(i,sortIndex);
            end%for
        end%function
    end%methods
end%classdef

%% References
%
% # M. Herzog, D. Schick, P. Gaal, R. Shayduk, C. von Korff Schmising & M.
% Bargheer (2011). _Analysis of ultrafast X-ray diffraction data in a 
% linear-chain model of the lattice dynamics_. Applied Physics A, 106(3), 
% 489–499. doi:10.1007/s00339-011-6719-z