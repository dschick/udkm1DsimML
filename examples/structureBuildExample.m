%% Structure Build Example
% This is a first example to get used to the basic element of each
% simulation of the udkm1dsim toolbox: a sample structure.
%
% For the 1D case a sample structure consists of a linear chain of
% different "layers". These layers themselves are called sub-structures and
% can again consist of sub-layers or finally of a certain number of unit
% cells. Each unit cell contains a number of atoms at certain relative
% positions in this unit cell.
%
% In this example script the basic concepts of creating atoms, unit cells,
% and sample structures are introduced. Furthermore one should easily see
% how to access all the properties of these physical objects.
%
% Be sure to include all folders of the udkm1dsim toolbox to your MATLAB 
% search path. Especially the _/parameters_ and _/helpers_ folders with all 
% their subfolders.

%%
% here we clear the workspace and close all figures
clear all;
close all;
% these are some units and constants we can use later on
u            = units;
const        = constants;

%% Atoms
% All necessary data are loaded from the paramter files defined by the given
% element symbol. As an additional 2nd parameter one can give an ID if you
% work with different atoms of the same element. The optinal 3rd input
% parameter is the ionicity of the atom, which has to be present in the
% parameter files.
O   = atomBase('O');
Op1 = atomBase('O','Op1',-1);
Op2 = atomBase('O','Op2',-2);
Ti  = atomBase('Ti');
Sr  = atomBase('Sr');
Ru  = atomBase('Ru');
Pb  = atomBase('Pb');
Zr  = atomBase('Zr');

% the available parameters of each atom can be easily displayed:
O.disp();
%%
Ru.disp();
%%
% Of course you can also access single properties, e.g. the mass of Ti atom:
disp(Ti.mass/u.kg);

%% Mixed Atoms
% If you want to create a mixed atom for a solid solution you can easily
% achieve this by the following lines of code.
% The input for the initialization of the mixed atom are the symbol, ID and
% name, whereas only the first is required.
ZT = atomMixed('ZT', 'ZT', 'Zircon-Titan 0.2 0.8');
ZT.addAtom(Zr, 0.2);
ZT.addAtom(Ti, 0.8);

% Let's have a look at the mixed properties:
ZT.disp();

%% Unit Cells
% Now we have created enough atoms to build a unit cell. We choose a cubic 
% unit cell of the perovskite SrRuO3 on a substrate of SrTiO3. 
% However we need to define some  properties of the unit cell first,
% Note that the thermal conductivity, linear thermal expansion, and heat
% capacity are defined either as constants or anonymous functions of the
% termperature T in Kelvin.

% c-axis lattice constants of the two materials
cSTOsub     = 3.905     *u.ang;
cSRO        = 3.94897   *u.ang;
% sound velocities [nm/ps] of the materials
svSRO       = 6.312     *u.nm/u.ps;
svSTO       = 7.800     *u.nm/u.ps;

%% SrRuO3 unit cell
propSRO.aAxis           = cSTOsub;              % aAxis
propSRO.bAxis           = cSTOsub;              % bAxis
propSRO.debWalFac       = 0;                    % Debye-Waller factor
propSRO.soundVel        = svSRO;                % sound velocity
propSRO.optPenDepth     = 43.8*u.nm;            % optical penetration depth
propSRO.thermCond       = 5.72*u.W/(u.m *u.K);  % heat conductivity
propSRO.linThermExp     = 1.03e-5;              % linear thermal expansion
propSRO.heatCapacity    = @(T)(455.2 + 0.112.*T - 2.1935e6./T.^2);
                                                % heat capacity [J/kg K]                              
% lets create the unit cell
SRO = unitCell('SRO','Strontium Ruthenate',cSRO,propSRO);
% add some atoms at relative positions in the unit cell
SRO.addAtom(O,0);
SRO.addAtom(Sr,0);
SRO.addAtom(O,0.5);
SRO.addAtom(O,0.5);
SRO.addAtom(Ru,0.5);

%% SrTiO3 substrate
propSTOsub.aAxis           = cSTOsub;           % aAxis
propSTOsub.bAxis           = cSTOsub;           % bAxis
propSTOsub.debWalFac       = 0;                 % Debye-Waller factor
propSTOsub.soundVel        = svSTO;             % sound velocity
propSTOsub.optPenDepth     = Inf;               % optical penetration depth
propSTOsub.thermCond       = 12*u.W/(u.m *u.K); % heat conductivity
propSTOsub.linThermExp     = 1e-5;              % linear thermal expansion
propSTOsub.heatCapacity    = @(T)(733.73 + 0.0248.*T - 6.531e6./T.^2);
                                                % heat capacity [J/kg K]

%% Non-Linear Strain Dependence
% In general the position of each atom in a unit cell dependence linearly
% from an external strain. In some cases this linear behaviour had to be 
% altered. This can be easily achieved by providing a strain dependent 
% anonymous function handle for the atom position when the atom is added to
% the unit cell.
STOsub = unitCell('STOsub', 'Strontium Titanate Substrate', cSTOsub, propSTOsub);
STOsub.addAtom(O, @(strain)(0.1*strain^2)); % quadratic strain dependency
STOsub.addAtom(Sr,0);
STOsub.addAtom(O,0.5);
STOsub.addAtom(O,0.5);
STOsub.addAtom(Ti,0.5);

% again we can display all properties of the unit cells:
STOsub.disp();
%%
% we can also visualize the positions of the atoms in the unit cell:
SRO.visualize();

%% Clone to Multiple
% For some simulations it may save computational time if a larger spacial
% grid is used. This can be achieved by e.g. providing larger unit cells.
% In order to automatically do that task we can easily clone unit cells to
% multiples of them:
SRO2 = SRO.clone2multiple(2);
SRO2.visualize();

%% Structure
% Since we have build two different kinds of unit cells we can build a
% sample structure
S = structure('Single Layer');
S.addSubStructure(SRO,100);     % add 100 layers of SRO to sample
S.addSubStructure(STOsub,1000); % add 1000 layers of STO substrate

% display the properties of the structure and visualize
S.disp();
%%
S.visualize();
%%
% there are various methods to get informations from the structure which
% are mostly self explanatory:
[dStart, dEnd, dMid] = S.getDistancesOfUnitCells(); 
K = S.getNumberOfSubStructures();
L = S.getNumberOfUniqueUnitCells();
M = S.getNumberOfUnitCells();
P = S.getAllPositionsPerUniqueUnitCell();
I = S.getDistancesOfInterfaces();
cAxis = S.getUnitCellPropertyVector('cAxis');

% we can also build more complicated structures like super lattices:
S2 = structure('Super Lattice');
% define a single double layer
DL = structure('Double Layer');
DL.addSubStructure(SRO,15);     % add 15 layers of SRO
DL.addSubStructure(STOsub,20);  % add 20 layers of STO substrate
% add the double layer to the super lattice
S2.addSubStructure(DL,10);      % add 10 double layers to super lattice
S2.addSubStructure(STOsub,500); % add 500 layers of STO substrate

% display and visualize the structure:
S2.disp();
%%
S2.visualize();
%%
% now there are some more interfaces than before:
I2 = S2.getDistancesOfInterfaces();
disp(length(I2));