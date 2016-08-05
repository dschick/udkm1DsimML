%% atomMixed
% The atomMixed class is sub class of atomBase and enables mixed atoms for
% certain alloys and stochiometric mixtures. All properties of the included
% sub-atoms of class atomBase are averaged and weighted with their
% stochiometric ratio
%
% Copyright (c) 2013, Daniel Schick, André Bojahr, Marc Herzog, Roman Shayduk, Clemens von Korff Schmising
% All rights reserved.
%
% License: BSD (use/copy/change/redistribute on own risk, mention the authors)

%% Classdef
% The atomMixture class inherited from atomBase
classdef atomMixed < atomBase
    %% Properties
    properties (SetAccess=public,GetAccess=public)
        atoms       = {};   % CELL ARRAY of atoms and fractions
        numAtoms    = 0;    % INTEGER number of atoms in AtomMixed
    end%properties
    
    %% Methods
    methods
        %% Constructor
        % Is executed each time an instance of this class is created. Only
        % the _symbol_ input is obligatory. 
        function obj = atomMixed(symbol,varargin)
            % initialize input parser and define defaults and validators
            p = inputParser;
            p.addRequired('symbol'        , @ischar);
            p.addOptional('ID'    , symbol, @ischar);
            p.addOptional('name'  , symbol, @ischar);
            % parse the input
            p.parse(symbol,varargin{:});
            % assign parser results to object properties
            % execute the super-class constructor here
            obj = obj@atomBase(p.Results.symbol,p.Results.ID);
            obj.name = p.Results.name;
        end%function
        
        %% addAtom
        % Add a atomBase instance with its stochiometric fraction to the
        % atomMixed instance
        function addAtom(obj,atom,fraction)
            obj.atoms(end+1,:) =  {atom fraction};
            obj.numAtoms = obj.numAtoms + 1;
            % calculate the mixed atomic properties of the atomMixed
            % instance
            obj.atomicNumberZ    = obj.atomicNumberZ     + fraction*   atom.atomicNumberZ;
            obj.massNumberA      = obj.massNumberA       + fraction*   atom.massNumberA;
            obj.mass             = obj.mass              + fraction*   atom.mass;
            obj.ionicity         = obj.ionicity          + fraction*   atom.ionicity;
        end%function
        
        %% Display
        % This method is called to display informations of the instance.
        function disp(obj)
            % display the properties of the atom
            disp@atomBase(obj);
            disp([num2str(obj.numAtoms) ' Constituents:']);
            for i = 1:obj.numAtoms
                fprintf('\t %s \t %3.2f%%\n', obj.atoms{i,1}.name, obj.atoms{i,2}*100);
            end%for
        end%function
        
        %% getAtomicFormFactor
        % Returns the mixed energy dependent atomic form factor.
        function f = getAtomicFormFactor(obj,E)
            f = 0;
            for j = 1:obj.numAtoms
                f = f + obj.atoms{j,1}.getAtomicFormFactor(E) * obj.atoms{j,2};
            end%for
        end%function
        
        %% getCMAtomicFormFactor
        % Returns the mixed energy and angle dependent atomic form factor.
        function f = getCMAtomicFormFactor(obj,E,qz)
            f = 0;
            for j = 1:obj.numAtoms
                f = f + obj.atoms{j,1}.getCMAtomicFormFactor(E,qz) * obj.atoms{j,2};
            end%for
        end%function
    end%methods
end%classdef