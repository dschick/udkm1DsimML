%% simulation
% The simulation class is the super class for all simulation classes.
%
% Copyright (c) 2013, Daniel Schick, André Bojahr, Marc Herzog, Roman Shayduk, Clemens von Korff Schmising
% All rights reserved.
%
% License: BSD (use/copy/change/redistribute on own risk, mention the authors)

%% Classdef
% Each simulation instance and all inherited class objects follow handle
% semantics. Hence a copy of such object will not copy the object itself,
% but only a handle to that object.
classdef simulation < handle
    %% Properties
    properties (SetAccess=public,GetAccess=public)
        S                           % OBJECT structure to simulate the phohon dynamics on
        forceRecalc     = false;    % BOOLEAN if true, everything is calculated despite of any saved old data
        cacheDir        = './';     % STRING path to cached data
        dispMessages    = true;     % BOOLEAN is true to display messages of from with in the simulations
        dispCalcTime    = true;     % BOOLEAN is true to display the duration of certain calculations 
                                    % (works only if displayMessages == true) 
        progressBarType = 'text';   % STRING type of the progressbar 'none', 'text', 'gui'
    end%properties
    %% Methods
    methods
        %% Constructor
        % Is executed each time an instance of this class is created. Only
        % the _structure_ input is obligatory.
        function obj = simulation(structure,forceRecalc,varargin)
            % initialize input parser and define defaults and validators
            p = inputParser;
            p.KeepUnmatched = true;
            p.addRequired('structure'           ,        @(x)isa(x,'structure'));
            p.addRequired('forceRecalc'         ,        @islogical);
            p.addParamValue('dispMessages'      , true,  @islogical);
            p.addParamValue('dispCalcTime'      , true,  @islogical);
            p.addParamValue('progressBarType'   ,'text', @(x)(ischar(x) & find(strcmp(x,{'none', 'text', 'gui'}))));
            p.parse(structure,forceRecalc,varargin{:});
            % assign parser results to object properties
            obj.S               = p.Results.structure;
            obj.forceRecalc     = p.Results.forceRecalc;
            obj.dispMessages    = p.Results.dispMessages;
            obj.dispCalcTime    = p.Results.dispCalcTime;
            obj.progressBarType = p.Results.progressBarType;
        end%function
        
        %% Display
        % This method is called to display informations of the instance.
        function disp(obj)
            disp('This is the current structure for the simulations:');
            disp('__________________________________________________');
            obj.S.disp();
            disp('__________________________________________________');
            disp('Display properties');
            disp(['force recalc             : ' bool2str(obj.forceRecalc)]);
            disp(['cache directory          : ' obj.cacheDir]);
            disp(['display messages         : ' bool2str(obj.dispMessages)]);
            disp(['display calculation time : ' bool2str(obj.dispCalcTime)]);
            disp(['progress bar type        : ' obj.progressBarType]);     
        end%function
                
        %% dispMessage
        % Displays the input message or the input message and input time.
        function dispMessage(obj,message,time)
            if nargin < 3 % no time is given, so its just a message
                if obj.dispMessages
                    disp(message);
                end%if
            else % this is a message with calculation time
                if obj.dispMessages && obj.dispCalcTime
                    disp([message ' ' num2str(time) ' seconds.']);
                end%if
            end%if
        end%function
        
        %% progressBar
        % Shows a progress bar depending on the value of
        % _obj.progressBarType_.
        function progressBar(obj,input)
            persistent h text; % remember the state
            if strcmp(obj.progressBarType,'text')
                % call external textprogressbar
                textprogressbar(input);
            elseif strcmp(obj.progressBarType,'gui')
                if isempty(h) && ~ischar(input),
                    % Progress bar must be initialized with a string
                    error('The text progress must be initialized with a string');
                elseif isempty(h) && ischar(input),
                    % Progress bar - initialization
                    h = waitbar(0,input); % open a waitbar
                    text = input;
                elseif ~isempty(h) && ischar(input),
                    % Progress bar  - termination
                    close(h);
                    clear text h;
                elseif ~isempty(h) && isnumeric(input)
                    % Progress bar - normal progress
                    waitbar(input/100,h, sprintf('%s %03.1f%%',text,input)); % it's just a waitbar
                else
                    % Any other unexpected input
                    error('Unsupported argument type');
                end%if
            end%if
        end%function
        
        %% odeProgressBar
        % Provides an interface for a progress bar for a ODE solver.
        % The progressbar can be 'text', 'gui' or 'none'.
        function status = odeProgressBar(obj,t,~,flag)
            persistent tf;
            
            if isempty(flag)
                % Integration steps
                ts=mean(t);
                progress=100*ts/tf;
                obj.progressBar(progress);
                status = 0;
            else
                switch flag
                    case 'init'     % Initializing progress bar
                        tf=max(t);
                        obj.progressBar('ODE integration: ');
                    case 'done'     % Finishing status function
                        tf=[];
                        obj.progressBar('');
                    otherwise
                        error('odetpbar:UnknownError',...
                            'Unknown error has occured');
                end%switch
            end%if
        end%function
        
        %% setCacheDir
        % Sets the path where calculation results can be saved for later
        % usage.
        function setCacheDir(obj, path)
            if exist(path,'dir')
                obj.cacheDir = path;
            else
                error('Path does not exist');
            end%if
        end%function 
    end%methods
end%classdef