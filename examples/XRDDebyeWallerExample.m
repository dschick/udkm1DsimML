%% XRD Example
% In this example we show how to carry out XRD simulations for the static
% and transient case.
%
% Before we can start a XRD simulation we need to build a sample
% structure. For the transient XRD simulations also according lattice 
% dynamics are necessary as input.
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

% create atoms
O   = atomBase('O');
Ti  = atomBase('Ti');
Sr  = atomBase('Sr');
Ru  = atomBase('Ru');
Pb  = atomBase('Pb');
Zr  = atomBase('Zr');

% c-axis lattice constants of the two materials
cSTOsub     = 3.905     *u.ang;
cSRO        = 3.94897   *u.ang;
% sound velocities [nm/ps] of the materials
svSRO       = 6.312     *u.nm/u.ps;
svSTO       = 7.800     *u.nm/u.ps;

%% SrRuO3 unit cell
propSRO.aAxis           = cSTOsub;              % aAxis
propSRO.bAxis           = cSTOsub;              % bAxis
propSRO.debWalFac       = {@(T)(0), @(T)(1e-10.*T)}; % Debye-Waller factor
propSRO.soundVel        = svSRO;                % sound velocity
propSRO.optPenDepth     = 43.8*u.nm;            % optical penetration depth
propSRO.thermCond       = {0,...                % electronic heat conductivity
                           0};                  % lattice heat conductivity [W/m K]
propSRO.linThermExp     = {1.03e-5, ...         % electronic linear thermal expansion
                           1.03e-5};            % lattice linear thermal expansion [1/K]
propSRO.heatCapacity    = {@(T)(0.112.*T), ...
                           @(T)(455.2 - 2.1935e6./T.^2)} ;
                                                % electronic and lattice
                                                % heat capacity [J/kg K]
propSRO.subSystemCoupling = {@(T)( 5e17*(T(2)-T(1)) ), ...
                             @(T)( 5e17*(T(1)-T(2)) )};
                                                % energy flux between
                                                % subsystems [W/m^3]




SRO = unitCell('SRO','SRO',cSRO,propSRO);
SRO.addAtom(O,0);
SRO.addAtom(Sr,0);
SRO.addAtom(O,0.5);
SRO.addAtom(O,0.5);
SRO.addAtom(Ru,0.5);


SRO.intLinThermExp    = {@(T)(1.03e-5.*T), @(T)(1.03e-5.*T)};
SRO.intHeatCapacity    = {@(T)(0.5*0.112.*T.^2), @(T)(455.2.*T + 0.112.*T.^2 - 2.1935e6./T)};

% SrTiO3 substrate
propSTOsub.aAxis           = cSTOsub;           % aAxis
propSTOsub.bAxis           = cSTOsub;           % bAxis
propSTOsub.debWalFac       = {@(T)(0), @(T)(1e-10.*T)};                 % Debye-Waller factor
propSTOsub.soundVel        = svSTO;             % sound velocity
propSTOsub.optPenDepth     = Inf;               % optical penetration depth
propSTOsub.thermCond       = {0,...             % electronic heat conductivity
                              0};               % lattice heat conductivity [W/m K]
propSTOsub.linThermExp     = {1e-5, ...         % electronic linear thermal expansion
                              1e-5};            % lattice linear thermal expansion [1/K]
propSTOsub.heatCapacity    = {@(T)(0.0248.*T), ...
                              @(T)(733.73 - 6.531e6./T.^2)} ;
                                                % electronic and lattice
                                                % heat capacity [J/kg K]
propSTOsub.subSystemCoupling = {@(T)( 5e17*(T(2)-T(1)) ), ...
                                @(T)( 5e17*(T(1)-T(2)) )};
                                                % energy flux between
                                                % subsystems [W/m^3]
                                         
STOsub = unitCell('STOsub', 'STOsub', cSTOsub, propSTOsub);
STOsub.addAtom(O,0);
STOsub.addAtom(Sr,0);
STOsub.addAtom(O,0.5);
STOsub.addAtom(O,0.5);
STOsub.addAtom(Ti,0.5);


STOsub.intLinThermExp    = {@(T)(1.03e-5.*T), @(T)(1.03e-5.*T)};
STOsub.intHeatCapacity    = {@(T)(0.5*0.112.*T.^2), @(T)(455.2.*T + 0.112.*T.^2 - 2.1935e6./T)};

% build the structure
S = structure('Single Layer');
S.addSubStructure(SRO,100);     % add 100 layers of SRO to sample
S.addSubStructure(STOsub,1000);  % add 250 layers of STO substrate

distances = S.getDistancesOfUnitCells();
%% Initialize Heat Simulation
cacheDir = './cache';
forceRecalc = false; 


%% set the parameters for the XRD simulations
forceRecalc     = false;
cacheDir        = './cache';
E               = 8047*u.eV; % X-ray energy
pol             = 0.5; % mixed X-ray polarization
theta           = (22:0.0001:24)*u.deg; % theta-range

% initialize the kinematical XRD simulation
K = XRDkin(S,forceRecalc,E,pol);
K.setCacheDir(cacheDir); % set the cache directory 
K.setQzByTheta(theta); % set the theta range
% calculate the static rocking curve for the homogenous sample in 1 line
Rkin = K.homogeneousReflectivity();

% initialize the dynamical XRD simulation
D = XRDdynDebyeWaller(S,forceRecalc,E,pol);
D.setCacheDir(cacheDir); % set the cache directory
D.setQzByTheta(theta); % set the theta range
% calculate the static rocking curve for the homogenous sample in 1 line
Rdyn = D.homogeneousReflectivity();
%% Static Rocking Curves - thin sample
% For the this sampple, both kinematical and dynamical XRD theory nealry
% coincident. However, the slight offset of the two dataset is due to the
% neglect of refraction in the kinematical model. The oscillations on both
% curves are so-called Kiessig fringes which originate from the finite
% thickness of the sample, which is smaller than the absorption/extinction
% depth of the X-ray photons.
figure(1);
semilogy(theta/u.deg,Rkin, '-k'); hold on;
semilogy(theta/u.deg,Rdyn, '-r');
ylim([1e-7 1]);
xlabel('\theta [deg.]');
ylabel('Reflectivity');
title('Static Rocking Curve');
legend('kinematical', 'dynamical');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');

%% Add static substrate
% In order to simulate rocking curves of samples having a quasi-infinite
% substrate, one can add such a static substrate to the sample structure as
% an special substructure. This substrate substructure is not considered in
% the _heat_ and _phonon_ simulations and only acts in XRD simulations.
substrate = structure('STOsubstrate');
substrate.addSubStructure(STOsub,100000);
% add the substrate to the sample
S.addSubstrate(substrate);

%% Recalculate the static rocking curves
Rkin = K.homogeneousReflectivity();
Rdyn = D.homogeneousReflectivity();

%% Static Rocking Curves - thick sample
% Obviously, the kinematical simulation fails at maximum of the substrate
% peak, since reflectivities greater than 1 are simulated. Furthermore, the
% kinematical result still shows Kiessig fringes which are vanished in the
% dynamical simulation. Both effects can be explained by the neglect of
% absorption in the kinematical approach. However, only dynamical theory is
% capable to calculate the Darwin width of the intensie subsrtate peak
% correctly, which is observable in the strong zoom of this Bragg peak.
figure(2);
subplot(2,1,1);
semilogy(theta/u.deg,Rkin, '-k'); hold on;
semilogy(theta/u.deg,Rdyn, '-r');
ylim([1e-7 10000]);
xlabel('\theta [deg.]');
ylabel('Reflectivity');
title('Static Rocking Curve incl. Static Substrate');
legend('kinematical', 'dynamical');
box on;
subplot(2,1,2);
semilogy(theta/u.deg,Rdyn, '-r');
ylim([1e-4 1]);
xlim([23.2 23.28]);
xlabel('\theta [deg.]');
ylabel('Reflectivity');
title('Darwin Width');
box on;
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');
%%
H = heat(S,forceRecalc);
H.setCacheDir(cacheDir); % set the cache directory
H.heatDiffusion = true;  % enable heat diffusion

time        = (-5:0.01:15)*u.ps;
fluence     = 15*u.mJ/u.cm^2;
pulseWidth  = 0.25*u.ps;
initTemp    = 300*u.K;

excitation(1,:) = fluence;    % fluence
excitation(2,:) = 0;          % time when the excitation happens
excitation(3,:) = pulseWidth; % pulse width of the excitation

% The resulting temperature profile is calculated in one line:
[tempMap deltaTempMap] = H.getTempMap(time,excitation,initTemp);

%% Plot the temperatures of each subsystem
figure(1)
subplot(2,1,1)
kk = surf(distances/u.nm,time/u.ps,tempMap(:,:,1));
set(kk, 'LineStyle', 'none');
title('Electron Temperature [K]'); 
xlabel('z [nm]'); ylabel('Delay [ps]');
axis([distances(1)/u.nm distances(end)/u.nm time(1)/u.ps time(end)/u.ps])
box; colorbar; colormap(jet(255));
caxis([300 3000]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');

subplot(2,1,2)
kk = surf(distances/u.nm,time/u.ps,tempMap(:,:,2));
set(kk, 'LineStyle', 'none');
title('Lattice Temperature [K]'); 
xlabel('z [nm]'); ylabel('Delay [ps]');
axis([distances(1)/u.nm distances(end)/u.nm time(1)/u.ps time(end)/u.ps])
box; colorbar; colormap(jet(255));
caxis([300 3000]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');


%% Plot the average electron and lattice temperature of the SRO layer
figure(2)
plot(time/u.ps, mean(tempMap(:,S.getAllPositionsPerUniqueUnitCell{1},1),2), '-r', 'LineWidth', 2);
hold on;
plot(time/u.ps, mean(tempMap(:,S.getAllPositionsPerUniqueUnitCell{1},2),2), '-k', 'LineWidth', 2);
title('Subsystem temperature');
xlabel('Delay [ps]');
ylabel('Temperature [K]');
box on; grid on;
ylim([0 3000]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');
legend('T_{electrons}', 'T_{lattice}');

%% Initialize Phonon Simulation
P = phononNum(S,forceRecalc);
P.setCacheDir(cacheDir); % set the cache directory
%% 
% Calculate coherent phonon dynamics (strainMap) providing the _time_, and 
% temperature inputs of the heat simulation:
strainMap = P.getStrainMap(time,tempMap,deltaTempMap);
% calculate a reduced number of strains per unique unit cell in order to
% save computational time for the transient XRD simulation
strainVectors = P.getReducedStrainsPerUniqueUnitCell(strainMap);
%%
% plot the results:
figure(4)
kk = surf(distances/u.nm,time/u.ps,strainMap);
set(kk,'LineStyle', 'none');
xlabel('z [nm]');
ylabel('Delay [ps]');
title('Strain - Coherent Phonon Dynamics');
axis([distances(1)/u.nm distances(end)/u.nm time(1)/u.ps time(end)/u.ps])
box on; colorbar; colormap(fireice(255));
caxis([-max(max(strainMap)) max(max(strainMap))]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');

%% Calculate transient XRD
% Reduce the number of angular steps, to save computational time
theta = (22:0.1:24)*u.deg;
D.setQzByTheta(theta);
% This is the parallel calculation of the transient XRD for the given 
% _strainMap_ 
R = D.getInhomogeneousReflectivity(strainMap,tempMap, 'sequential');

%%
% plot the results:
figure(5)
kk = surf(theta/u.deg,time/u.ps,log10(R));
set(kk,'LineStyle', 'none');
colorbar; colormap(fireice(255));
axis xy;
xlabel('\theta [deg.]');
ylabel('Delay [ps]');
title('Transient X-ray Reflectivity');
view(0,90); box on;
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');

%%
% convolute the results with an angular instrument function (Pseudo Voigt)
% and temporal instrument function (Gauss)
instFunc = @(theta)(pseudo_voigt(theta,0.02*u.deg));
Rconv = D.convWithInstrumentFunction(R,theta,instFunc);
Rconv = D.convWithTemporalResolution(Rconv,time,1*u.ps);

%%
% plot the convoluted results:
figure(6)
kk = surf(theta/u.deg,time/u.ps,log10(Rconv));
set(kk,'LineStyle', 'none');
colorbar; colormap(fireice);
axis xy;
xlabel('\theta [deg.]');
ylabel('Delay [ps]');
title('Transient X-ray Reflectivity convoluted in \theta and Delay');
view(0,90);
box on;
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');