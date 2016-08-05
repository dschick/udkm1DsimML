%% Phonon Anharmonic Example
% In this example we show how to include anharmonicity and damping for the
% calculation of the coherent phonon dynamics using the numerical phonon
% simulation model.
%
% Before we can start a phonon simulation we need to build a sample
% structure and simulate the temperature map after the excitation.
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
%% Build a Sample Structure
O   = atomBase('O');
Ti  = atomBase('Ti');
Sr  = atomBase('Sr');
Ru  = atomBase('Ru');
Pb  = atomBase('Pb');
Zr  = atomBase('Zr');

cSTOsub     = 3.905     *u.ang;
cSRO        = 3.94897   *u.ang;
svSRO       = 6.312     *u.nm/u.ps;
svSTO       = 7.800     *u.nm/u.ps;

propSRO.aAxis           = cSTOsub;              % aAxis
propSRO.bAxis           = cSTOsub;              % bAxis
propSRO.debWalFac       = 0;                    % Debye-Waller factor
propSRO.soundVel        = svSRO;                % sound velocity
propSRO.optPenDepth     = 43.8*u.nm;            % optical penetration depth
propSRO.thermCond       = 5.72*u.W/(u.m *u.K);  % heat conductivity
propSRO.linThermExp     = 1.03e-5;              % linear thermal expansion
propSRO.heatCapacity    = @(T)(455.2 + 0.112.*T - 2.1935e6./T.^2);
                                                % heat capacity [J/kg K]
SRO = unitCell('SRO','SRO',cSRO,propSRO);
SRO.addAtom(O,0);
SRO.addAtom(Sr,0);
SRO.addAtom(O,0.5);
SRO.addAtom(O,0.5);
SRO.addAtom(Ru,0.5);

propSTOsub.aAxis           = cSTOsub;           % aAxis
propSTOsub.bAxis           = cSTOsub;           % bAxis
propSTOsub.debWalFac       = 0;                 % Debye-Waller factor
propSTOsub.soundVel        = svSTO;             % sound velocity
propSTOsub.optPenDepth     = Inf;               % optical penetration depth
propSTOsub.thermCond       = 12*u.W/(u.m *u.K); % heat conductivity
propSTOsub.linThermExp     = 1e-5;              % linear thermal expansion
propSTOsub.heatCapacity    = @(T)(733.73 + 0.0248.*T - 6.531e6./T.^2);
                                                % heat capacity [J/kg K]
STOsub = unitCell('STOsub', 'STOsub', cSTOsub, propSTOsub);
STOsub.addAtom(O,0);
STOsub.addAtom(Sr,0);
STOsub.addAtom(O,0.5);
STOsub.addAtom(O,0.5);
STOsub.addAtom(Ti,0.5);
%% Anharmonicity & Damping
% We can easily change the phonon propagation in the numerical model by
% introducing phonon damping and higher order spring constants to any

STOsub.phononDamping = 1e-12;        % [kg/s]
STOsub.setHOspringConstants([-7e12]);% qubic potential [kg/m s^2]

S = structure('Single Layer');
S.addSubStructure(SRO,100);     % add 100 layers of SRO to sample
S.addSubStructure(STOsub,2000); % add 2000 layers of STO substrate

distances = S.getDistancesOfUnitCells();
%% Initialize Heat Simulation
cacheDir = './cache';
forceRecalc = false; 

H = heat(S,forceRecalc);
H.setCacheDir(cacheDir); % set the cache directory
time        = (-10:0.1:90)*u.ps;
fluence     = 30*u.mJ/u.cm^2;
initTemp    = 300*u.K;
%% Excitation
% Calculate the temperature map after instantaneous excitation at $t=0$.
% For the phonon simulations also the differential temperature map 
% _deltaTempMap_in time is necessary:
[tempMap deltaTempMap] = H.getTempMap(time,fluence,initTemp);
%% Numerical Phonon Simulation
clear P;
P = phononNum(S,forceRecalc);
P.setCacheDir(cacheDir); % set the cache directory
%% 
% The coherent phonon dynamics (strainMap) is calculated by one line
% providing the _time_, and temperature inputs of the heat simulation:
strainMap = P.getStrainMap(time,tempMap,deltaTempMap);
%%
% plot the results:
figure(1)
kk = surf(distances/u.nm,time/u.ps,strainMap);
set(kk, 'LineStyle', 'none');
title('Strain'); xlabel('z [nm]'); ylabel('Time [ps]');
axis([distances(1)/u.nm distances(end)/u.nm time(1)/u.ps time(end)/u.ps])
box; colorbar; colormap(fireice(255));
caxis([-max(max(strainMap)) max(max(strainMap))]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');
%%
figure(2)
show = 101:150:length(time);
plot(distances/u.nm,strainMap(show,:));
xlim([distances(1)/u.nm distances(end)/u.nm]);
xlabel('z [nm]'); ylabel('Strain');
legend([num2str(time(show)'/u.ps) repmat(' ps',length(show),1) ],'Location','EastOutside');
box on; grid on;
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');