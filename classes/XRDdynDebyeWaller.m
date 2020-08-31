%% XRDdyn
% The XRDdyn class simulates dynamical X-ray diffraction on a 1D structure.
%
% Copyright (c) 2013, Daniel Schick, André Bojahr, Marc Herzog, Roman Shayduk, Clemens von Korff Schmising
% All rights reserved.
%
% License: BSD (use/copy/change/redistribute on own risk, mention the authors)

%% Classdef
% Each XRDdynDebyeWaller instance and all inherited class objects are inherted from 
% the XRDdyn class which follows handle semantics. 
% Hence a copy of such object will not copy the object itself,
% but only a handle to that object.
classdef XRDdynDebyeWaller < XRDdyn
    %% Methods
    methods
        %% Constructor
        % Is executed each time an instance of this class is created. Only
        % the _structure_ input is obligatory.
        function obj = XRDdynDebyeWaller(structure,varargin)
            obj = obj@XRDdyn(structure,varargin{:});
        end%functions
        
        %% Display
        % This method is called to display informations of the instance.
        function disp(obj)
            disp('Dynamical X-Ray Diffraction Debye Waller simulation properties:');
            % call the parent display method
            disp@XRD(obj);
        end%function  
        
        %% getHash
        % Returns a unique hash given by the energy $E$, $q_z$ range,
        % polarization factor and the strain vectors as well as the sample
        % structure hash.
        function hash = getHash(obj,strainMap,tempMap)
            % dataHash is an external function
            if nargin == 2
                hash = [obj.S.getHash('XRD') '_' dataHash({obj.E obj.qz obj.pol})];
            else
                % reduce size of strainMap when it has more than 1e6 elements
                if numel(strainMap) > 1e6
                    strainMap = reshape(strainMap,1,numel(strainMap));
                    strainMap = strainMap(1:1e6);
                end
                if numel(tempMap) > 1e6
                    tempMap = reshape(tempMap,1,numel(tempMap));
                    tempMap = tempMap(1:1e6);
                end
                hash = [obj.S.getHash('XRD') '_' dataHash({obj.E obj.qz obj.pol strainMap tempMap})];
            end
        end%functions    
                
        %% getInhomogeneousReflectivity
        % Returns the reflectivity of an inhomogenously strained sample
        % structure for a given _strainMap_ in position and time, as well 
        % as for a given set of possible strains for each unit cell in the
        % sample structure (_strainVectors_).
        % If no reflectivity is saved in the cache it is caluclated.
        % Providing the _type_ (parallel [default], sequential,
        % distributed) for the calculation the corresponding subroutines
        % for the reflectivity computation are called:
        %
        % * *parallel* parallelization over the time steps utilizing 
        % MATLAB's Parallel Computing Toolbox
        % * *distributed* parallelization over the time steps utilizing 
        % MATLAB's Distribted Computing Toolbox
        % * *sequential* no parallelization at all
        function R = getInhomogeneousReflectivity(obj,strainMap,tempMap,varargin)
            % create a hash of all simulation parameters
            filename = ['inhomogeneousReflectivityDyn_' obj.getHash(strainMap, tempMap) '.mat'];
            fullfilename = fullfile(obj.cacheDir, filename);
            % check if we find some corresponding data in the cache dir
            if exist(fullfilename,'file') && ~obj.forceRecalc
                % found something so load it
                load(fullfilename);
                obj.dispMessage(['_inhomogeneousReflectivity_ loaded from file ' fullfilename]);
            else
                tic
                obj.dispMessage('Calculating _inhomogenousReflectivity_ ...');
                % parse the input arguments
                p = inputParser;
                p.KeepUnmatched = true;
                p.addRequired('strainMap'         , @isnumeric);
                p.addRequired('tempMap'         , @isnumeric);
                p.addOptional('type',   'parallel', @(x)(ischar(x) & find(strcmp(x,{'parallel', 'sequential', 'distributed'}))));
                p.addOptional('job' ,   '');
                p.addOptional('numWorker',   1    , @isnumeric);
                p.parse(strainMap,tempMap,varargin{:});
                % assign parser results to object properties
                type          = p.Results.type;
                strainMap     = p.Results.strainMap;
                tempMap     = p.Results.tempMap;
                job           = p.Results.job;
                numWorker     = p.Results.numWorker;
                N = obj.S.getNumberOfUnitCells; % nb of unit cells
                M = size(tempMap, 1);               % nb of time steps
                K = obj.S.numSubSystems;        % nb of subsystems
                tempMap = reshape(tempMap, [M, N, K]);
                
                % select the type of computation
                switch type
                    case 'parallel'
                        R = parallelInhomogeneousReflectivity(obj,strainMap,tempMap);
                    case 'distributed'
                        R = distributedInhomogeneousReflectivity(obj,strainMap,tempMap,job,numWorker);
                    otherwise % sequential
                        R = sequentialInhomogeneousReflectivity(obj,strainMap,tempMap);
                end%switch
                
                obj.dispMessage('Elapsed time for _inhomogeneousReflectivity_:',toc);
                save(fullfilename, 'R');
                obj.dispMessage(['_inhomogeneousReflectivity_ saved to file ' fullfilename]);
            end%if
        end%function
            
        %% parallelInhomogeneousReflectivity
        % Returns the reflectivity of an inhomogenously strained sample
        % structure for a given _strainMap_ in position and time, as well 
        % as for a given set of possible strains for each unit cell in the
        % sample structure (_strainVectors_).
        % The function tries to parallize the calculation over the time
        % steps (_parallel = true_, since the results do not depent on each 
        % other. The routine checks whether the MATLAB pool is open - if 
        % not it opens the matlab pool with the default configuration.
        function R = parallelInhomogeneousReflectivity(obj,strainMap,tempMap)
            %initialize
            N = size(strainMap,1); % time steps
            R = zeros(N,length(obj.qz));            
            
            if verLessThan('matlab', '8.5') % this is everything before MATLAB 2015a
                s = matlabpool('size'); % get the size of the matlabpool
                if s == 0 % no matlabpool open
                    obj.dispMessage(['No matlab pool was opened in advance, so lets do it now with the default configuration!']);
                    matlabpool open;
                end%if
            else % this is for everthing starting with MATLAB 2015a
                if isempty(gcp('nocreate')) %s == 0 % no matlabpool open
                    obj.dispMessage(['No matlab pool was opened in advance, so lets do it now with the default configuration!']);
                    parpool;
                end%if
            end%if

            % check for path of ParforProgMon class to add it to
            % javapath
            p = fileparts(which('ParforProgMon.m'));
            str = ['javaaddpath ' p];
            feval(@pctRunOnAll,str);

            % make progresspar with the external parforProgressMonitor
            % package
            ppm = ParforProgMon('Please wait... ',N);
            parfor i = 1:N
                ppm.increment();
                % get the inhomogenous reflectivity of the sample
                % structure for each time step of the strain map
                R(i,:) = obj.calcInhomogeneousReflectivity(strainMap(i,:),tempMap(i,:,:));
            end%parfor
        end%function
        
        %% sequentialInhomogeneousReflectivity
        % Returns the reflectivity of an inhomogenously strained sample
        % structure for a given _strainMap_ in position and time, as well 
        % as for a given set of possible strains for each unit cell in the
        % sample structure (_strainVectors_).
        % The function calculates the results sequentially without
        % parallelization.
        function R = sequentialInhomogeneousReflectivity(obj,strainMap,tempMap)
            %initialize
            N = size(strainMap,1); % time steps
            R = zeros(N,length(obj.qz));
            
            obj.progressBar('Please wait... '); % open a progress bar
            for i = 1:N
                % get the inhomogenous reflectivity of the sample
                % structure for each time step of the strain map
                R(i,:) = obj.calcInhomogeneousReflectivity(strainMap(i,:),tempMap(i,:,:));
                % print the progress to console
                obj.progressBar(i/N*100);
            end%for
            obj.progressBar('');                
        end%function
        
        %% distributedInhomogeneousReflectivity
        % Return the reflectivity of an inhomogenous sample structure for 
        % a given _strainMap_ in position and time, as well as for a given 
        % set of possible strains for each unit cell in the sample 
        % structure (_strainVectors_). This method is distributed over 
        % several workers using the MATLAB dist. computing toolbox. It 
        % requires a Job handle and the number of workers that should 
        % contribute.                
        function R = distributedInhomogeneousReflectivity(obj,Job,numWorker,strainMap,tempMap)
            % initialize
            N           = size(strainMap,1);
            taskSize    = floor(N/numWorker);
            rest        = mod(N,numWorker);
            numTasks    = 0;                

            % traverse all tasks
            for i = 1:taskSize:(taskSize*numWorker)
                % create a task for each part of the strain pattern
                createTask(Job, @obj.getInhomogeneousReflectivity, 1, {strainMap(i:(i+taskSize-1),:) tempMap(i:(i+taskSize-1),:,:) 'sequential'});
                numTasks = numTasks+1;
            end%for
            % if there are parts left in the strain pattern, we have to
            % distribute them no at the end
            if rest > 0
                for i = taskSize*numWorker:taskSize:(taskSize*numWorker+rest)
                    if i+taskSize > N
                        i_end = N;
                    else
                        i_end = i+taskSize;
                    end%if
                    % create a task for each part of the left strain pattern
                    createTask(Job, @obj.parallelInhomogeneousReflectivity, 1, {strainMap(i+1:i_end,:) tempMap(i+1:i_end,:,:) 0});
                    numTasks = numTasks+1;
                end%for
            end%if

            %Run the job.
            submit(Job);
            obj.dispMessage('Job submitted. Waiting for tasks ...');
            % plot the list of tasks
            get(Job, 'Tasks')

            % Wait for the job to finish running tasks, and retrieve the job results.
            obj.progressBar('Please wait... ');
            while ~waitForState(Job,'finished',1)
                [~, ~, finished] = findTask(Job);
                obj.progressBar(length(finished)/numTasks*100);
            end%while
            obj.progressBar('');

            % get the results
            out = getAllOutputArguments(Job);
            % build the output reflectifity from the distributed results
            R = zeros(N,length(obj.qz));
            m = 1;
            for i = 1:taskSize:(taskSize*numWorker)
                R(i:(i+taskSize-1),:) = out{m};
                m = m+1;
            end%for
            % if we have a rest, we add it now
            if rest > 0
                for i = taskSize*numWorker:taskSize:(taskSize*numWorker+rest)
                    if i+taskSize > N
                        i_end = N;
                    else
                        i_end = i+taskSize;
                    end
                    R(i+1:i_end,:) = out{m};
                    m = m+1;
                end%for
            end%if
        end%function
         
        %% calcInhomogeneousReflectivity
        % Calculates the reflectivity of a inhomogenous sample structure 
        % for a given strain vector for a single time step. Similar to the
        % homogeneous sample structure, the reflectivity of an unit cell is
        % calculated from the reflection-transmission matrices $H_i$ of 
        % each atom and the phase matrices between the atoms $L_i$:
        %
        % $$ M_{RT} = \prod_i H_i \, L_i $$
        %
        % Since all layers are generally inhomogeneously strained we have 
        % to traverse all individual unit cells ($j = 1\ldots M$) in the 
        % sample to calculate the total reflection-transmission matrix 
        % $M_{RT}^t$:
        %
        % $$ M_{RT}^t = \prod_{j=1}^M M_{RT,j} $$
        % 
        % The reflectivity of the $2\times 2$ matrices for 
        % each $q_z$ is calculates as follow:
        %
        % $$ R = \left|M_{RT}^t(1,2)/M_{RT}^t(2,2)\right|^2 $$
        %
        function R = calcInhomogeneousReflectivity(obj,strains,temps)            
            % if no all-ref-trans matrices are given, we have to calculate
            % them first.
            % initialize
            N = obj.S.getNumberOfUnitCells();
            RT = repmat(eye(2,2),[1 1 length(obj.qz)]);
            % traverse all unitCells in the sample structure
            for i = 1:N
                % Find the ref-trans matrix in the RTM cell array for the
                % current unitCell ID and applied strain. Use the
                % _knnsearch_ funtion to find the nearest strain value.
                uc = obj.S.getUnitCellHandle(i);
                RT = mtimesx(RT,obj.getUCRefTransMatrix(uc,strains(i),squeeze(temps(:,i,:))));
            end%for
            % add the reflectivity of a substrate of available
            if ~isempty(obj.S.substrate)
                RT = mtimesx(RT,obj.homogeneousRefTransMatrix(obj.S.substrate));
            end
            % calculate reflectivity from ref-trans matrix
            R = obj.getReflectivityFromMatrix(RT);
        end%function
    end%methods
end%classdef

%% References
%
% # J. Als-Nielson, & D. McMorrow (2001). _Elements of Modern X-Ray 
% Physics_. New York: John Wiley & Sons, Ltd. doi:10.1002/9781119998365