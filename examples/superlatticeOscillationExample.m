%% Superlattice Oscillations
% In this example a superlattice (SL) structure of 11 DL consisting of 20 
% SRO layers and 38 STO layers on a STO substrate are excited by a single
% ultrashort laser pulse. The excitation of an optical SL phonon mode
% results in the intensity oscillation of according Bragg peaks.
%
% The simulation can be compared to the experiments given in Ref. [1].

%% general setting
close all; % close all previous figures
clear all; % clear the previous workspace
cacheDir = './cache'; % set the cache folder
forceRecalc = false;   % force the simulations to recalc

% define units and constants objects
u = units;
const = constants;

%% Set simulation parameters
time            = (-1:0.01:2)*u.ps;        % the time we want to simulate
E               = 500*u.eV;               % energy of X-rays
sp              = 0;                        % polarization factor: 
                                            % 0 -> S; 0.5 -> mixed; 1 -> P
qz              = (0.0:0.05:10)*u.ang^-1;% qz range
theta           = (0:0.01:45)*u.deg;
initTemp        = 300*u.K;                  % initial temperature of the sample
fluence         = 20*u.mJ/u.cm^2;
heatDiffusion   = false;                    % disable heat diffusion

%% Initialize the Sample
disp('Initialize sample ...');

%% Initialize the atoms of the Sample Structure
O  = atomBase('O');
Ti = atomBase('Ti');
Sr = atomBase('Sr');
Ru = atomBase('Ru');

%% Initialize the unitCells of the Sample
% lattice constants in Angstrom
cSTOsub     = 3.905     *u.ang;
cSTO        = 3.921     *u.ang;
cSRO        = 3.9493    *u.ang;

% sound velocities [nm/ps]
svSRO       = 6.312     *u.nm/u.ps;
svSTO       = 7.900     *u.nm/u.ps;

%%% SRO layer
propSRO.aAxis           = cSTOsub;              % aAxis
propSRO.bAxis           = cSTOsub;              % bAxis
propSRO.debWalFac       = 0;                    % Debye-Waller factor
propSRO.soundVel        = svSRO;                % sound velocity
propSRO.optPenDepth     = 52*u.nm;              % optical penetration depth [nm]
propSRO.thermCond       = 5.72*u.W/(u.m *u.K);  % heat conductivity [W/m K]
propSRO.linThermExp     = 1.03e-5;              % linear thermal expansion
propSRO.heatCapacity    = 464.419;              % heat capacity [J /kg K]

SRO = unitCell('SRO','SRO',cSRO,propSRO);
SRO.addAtom(O,0);
SRO.addAtom(Sr,0);
SRO.addAtom(O,0.5);
SRO.addAtom(O,0.5);
SRO.addAtom(Ru,0.5);

%%% STO layer
propSTO.aAxis           = cSTOsub;                  % aAxis
propSTO.bAxis           = cSTOsub;                  % bAxis
propSTO.debWalFac       = 0;                        % Debye-Waller factor
propSTO.soundVel        = svSTO;                    % sound velocity
propSTO.optPenDepth     = Inf;                      % optical penetration depth
propSTO.thermCond       = 12*u.W/(u.m *u.K);        % heat conductivity
propSTO.linThermExp     = 1e-5;                     % linear thermal expansion
propSTO.heatCapacity    = @(T)(733.73 + 0.0248.*T - 6.531e6./T.^2);
                                                % heat capacity [J/kg K]

STO = unitCell('STO', 'STO', cSTO, propSTO);
STO.addAtom(O,0);
STO.addAtom(Sr,0);
STO.addAtom(O,0.5);
STO.addAtom(O,0.5);
STO.addAtom(Ti,0.5);

%%% STO substrate
propSTOsub.aAxis           = cSTOsub;           % aAxis
propSTOsub.bAxis           = cSTOsub;           % bAxis
propSTOsub.debWalFac       = 0;                 % Debye-Waller factor
propSTOsub.soundVel        = svSTO;             % sound velocity
propSTOsub.optPenDepth     = Inf;               % optical penetration depth [nm]
propSTOsub.thermCond       = 12*u.W/(u.m *u.K); % heat conductivity [W/m K]
propSTOsub.linThermExp     = 1e-5;              % linear thermal expansion
propSTOsub.heatCapacity    = @(T)(733.73 + 0.0248.*T - 6.531e6./T.^2);
                                                % heat capacity [J/kg K]
STOsub = unitCell('STOsub', 'STOsub', cSTOsub, propSTOsub);
STOsub.addAtom(O,0);
STOsub.addAtom(Sr,0);
STOsub.addAtom(O,0.5);
STOsub.addAtom(O,0.5);
STOsub.addAtom(Ti,0.5);

%% Initialize the Sample Structure
DL = structure('DL: 20xSRO + 38xSTO');
DL.addSubStructure(SRO,12);
DL.addSubStructure(STO,11);

S = structure('11xDL on STO');
% add DL to the structure
S.addSubStructure(DL,100);
% S.addSubStructure(STOsub,1000);
% add a static substrate to the sample structure
substrate = structure('STOsubstrate');
substrate.addSubStructure(STOsub,1000000)% 
% S.addSubstrate(substrate);

distances = S.getDistancesOfUnitCells(); % these are the distances of each unitCell from the surface

%% Initialize the Heat Diffusion Simulation
% Provide the sample structure for the heat simulation.
H = heat(S,forceRecalc,heatDiffusion);
H.setCacheDir(cacheDir); % set the cache directory

%% Calculate Temperature Pattern
% Here calculate the temperature map and temperature difference map for the
% given time, fluence and initial temperature.
[tempMap, deltaTempMap] = H.getTempMap(time,fluence,initTemp);
%% Plot the Results of the Heat Simulations
figure(1)
% plot the temperature map
kk = surf(distances/u.nm,time/u.ps,tempMap);
set(kk, 'LineStyle', 'none');
xlabel('z [nm]'); 
ylabel('Delay [ps]');
title('Temperature [K]');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');
axis([distances(1)/u.nm distances(end)/u.nm time(1)/u.ps time(end)/u.ps]);
box on; colorbar('Location', 'EastOutside'); colormap jet(255);

%% Initialize the Phonon Dynamics Simulation
% Provide the sample structure for the phonon dynamics simulation
P = phononNum(S,forceRecalc);
P.setCacheDir(cacheDir); % set the cache directory

%% Calculate Strain Pattern
% Here we calculate the strain map for the given temperature, and 
% temperature difference map and the time vector.
strainMap = P.getStrainMap(time,tempMap,deltaTempMap);

%% Plot the Results of the Phonon Dynamics Simulations
figure(2);
kk = surf(distances/u.nm,time/u.ps,100*strainMap);
set(kk, 'FaceColor', 'interp', 'LineStyle', 'none');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');
axis([distances(1)/u.nm distances(end)/u.nm time(1)/u.ps time(end)/u.ps])
caxis([-100*max(max(strainMap)) 100*max(max(strainMap))]);
colormap fireice(255);
colorbar('location', 'EastOutside');
xlabel('z [nm]');
ylabel('Delay [ps]');
title('Strain [%]');

%% Initilaize dynamic XRD Simulation
disp('Initilaize dynamic XRD Simulation');
D = XRDdyn(S,forceRecalc,E,sp); % set main parameters
D.setCacheDir(cacheDir); % set the cache directory
% D.setQz(qz); % set q_z range
D.setQzByTheta(theta);

%% Homogeneous Dynamic XRD Simulations
% Here we calculate the rocking curve for no strain at all for the full 
% sample .
[Rh, A] = D.homogeneousReflectivity(); % thats all

% define an instrumental function to convolute the result with:
mu          = 0.1;
width       = 0.001*u.ang^-1;
instFunc    = @(qz)(pseudo_voigt(qz,width,mu));

% carry out the convolution with the instrumental function 
% [Rhi, xc] = D.convWithInstrumentFunction(Rh,qz,instFunc);

% Plot the Results of the homogeneous Dynamic XRD Simulations
figure(3)
semilogy(theta/u.deg,D.getReflectivityFromMatrix(A{1}{1})*10, '-r', 'LineWidth', 1);
hold on;
semilogy(theta/u.deg,D.getReflectivityFromMatrix(A{1}{2})*10, '-b', 'LineWidth', 1); 
semilogy(theta/u.deg,D.getReflectivityFromMatrix(A{2})*100, '-', 'Color',[0 0.5 0], 'LineWidth', 1);
semilogy(theta/u.deg, D.getReflectivityFromMatrix(A{3}), '-k', 'LineWidth', 2);
hold off
%axis([qz(1)/u.ang^-1 qz(end)/u.ang^-1 1e-5 1]);
set(gca,'YScale', 'log');
grid off; box on;
title('Homogeneous XRD');
xlabel('theta [deg]'); 
ylabel('Reflectivity');
lh = legend('20xSRO','38xSTO', '1xDL', '11xDL + Substrate');
set(lh, 'Location', 'NorthWest', 'FontSize', 10);
set(gca, 'Ydir', 'normal', 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');

%% Inhomogeneous Dynamic XRD Simulations
% Now we calculate the inhomogeneous dynamic XRD simulation for each time
% step each single strain for each unit cell in the sample

%%% 
% Calculate all strains that should be applied to all 
% unique types of unit cells in the sample structure in order to save
% computational time.
strainVectors = P.getReducedStrainsPerUniqueUnitCell(strainMap,100);

R = D.getInhomogeneousReflectivity(strainMap,strainVectors); % thats all

%% Plot the Results of the Dynamic XRD Simulations

% convolute with instrument function and temporal resolution first
%RI = D.convWithInstrumentFunction(R,qz,instFunc);
%RI = D.convWithTemporalResolution(RI,time,0.2*u.ps);

RI = R;

figure(4)
kk = surf(theta/u.deg,time/u.ps,log10(RI));
set(kk, 'LineStyle', 'none');
% axis([qz(1)/u.ang^-1 qz(end)/u.ang^-1 time(1)/u.ps time(end)/u.ps])
axis([theta(1)/u.deg theta(end)/u.deg time(1)/u.ps time(end)/u.ps])
box on; colorbar; colormap(fireice(256));
caxis([-6 0]);
xlabel('q_z [Ang^{-1}]'); 
ylabel('Delay [ps]');
title('Transient X-ray Reflectivity');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');

%% Plot transient Rocking Curve
figure(5)
semilogy(theta/u.deg,RI(time == 0.55*u.ps, :), '-r', 'LineWidth', 2);hold on;
semilogy(theta/u.deg,RI(time == -1*u.ps, :), '-k', 'LineWidth', 2);
axis([theta(1)/u.deg theta(end)/u.deg 1e-5 1]);
grid off; box on;
xlabel('q_z [Ang^{-1}]');
ylabel('Reflectivity');
title('Rocking Curves');
legend('t < 0', 't = 1.5 ps');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');
hold off;

%% Plot Bragg Peak Intensity Oscillations
% SL2 = sum(RI(:,qz > 3.195*u.ang^-1 & qz < 3.2*u.ang^-1),2);
SL0 = sum(RI(:,theta > 7.6*u.deg & theta < 9.4*u.deg),2);

figure(6)
subplot(2,1,1);
% plot(time/u.ps,SL2/SL2(time == -1*u.ps), '-b', 'LineWidth', 2);
% axis([-2 22 0 1.1]);
% title('Intensity SL 0');
% ylabel('\Delta R / R_0');
% set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');

subplot(2,1,2);
plot(time/u.ps,SL0/SL0(time == -1*u.ps), '-b', 'LineWidth', 2);
% axis([-2 22 0 50]);
title('Intensity SL +2');
xlabel('Delay [ps]');
ylabel('\Delta R / R_0');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');

%%%
disp('finished');

%% References
%
% # M. Herzog, D. Schick, W. Leitenberger, R. Shayduk, R. M. van der Veen, 
% C. J. Milne, S. L. Johnson, et al. (2012). _Tailoring interference and 
% nonlinear manipulation of femtosecond x-rays_. New Journal of Physics, 
% 14(1), 13004. doi:10.1088/1367-2630/14/1/013004