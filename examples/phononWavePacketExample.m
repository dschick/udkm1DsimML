%% Phonon Wave Packet
% In this example a 15nm metallic SRO thin film on dielectric STO substrate
% is excited by a multipulse sequence of ultrashort laser pulses with a 
% pulse distance of 7.2ps in order to generate a quasi-monochromatic phonon 
% wave packet. The phonon wave packet is evidenced as separate side-bands of the 
% substrate peak Bragg reflection.
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
time            = (-5:0.5:100)*u.ps;          % the time we want to simulate
E               = 12000*u.eV;               % energy of X-rays
sp              = 0;                        % polarization factor: 
                                            % 0 -> S; 0.5 -> mixed; 1 -> P
qz              = (3.1:0.0005:3.3)*u.ang^-1; % qz range
initTemp        = 300*u.K;                  % initial temperature of the sample
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

STOsub.phononDamping = 1e-12; % phonon damping [kg/s]

%% Initialize the Sample Structure
S = structure('SRO on STO');
% add unit cells to the structure
S.addSubStructure(SRO,38);
S.addSubStructure(STOsub,5000);
% add a static substrate to the sample structure
substrate = structure('STOsubstrate');
substrate.addSubStructure(STOsub,1000000)% 
S.addSubstrate(substrate);

distances = S.getDistancesOfUnitCells(); % these are the distances of each unitCell from the surface

%% Initialize the Heat Diffusion Simulation
% Provide the sample structure for the heat simulation.
H = heat(S,forceRecalc,heatDiffusion);
H.setCacheDir(cacheDir); % set the cache directory

%% Define the Multi-Pulse Excitation
pulseWidth      = 0*u.ps; 
pulseSep        = 7.2*u.ps; 
fluence         = 5.5*u.mJ/u.cm^2;

excitation(1,:) = fluence*ones(8,1); % fluence
excitation(2,:) = (0:7)*pulseSep; % time when the excitation happens
excitation(3,:) = pulseWidth;   % pulse width of the excitation

%% Calculate Temperature Pattern
% Here calculate the temperature map and temperature difference map for the
% given time, fluence and initial temperature.
[tempMap deltaTempMap] = H.getTempMap(time,excitation,initTemp);
%% Plot the Results of the Heat Simulations
figure(1)
% plot the temperature map
kk = surf(distances/u.nm,time/u.ps,tempMap);
set(kk, 'LineStyle', 'none');
xlabel('z [nm]'); 
ylabel('Delay [ps]');
title('Temperature [K]');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');
axis([0 100 time(1)/u.ps time(end)/u.ps]);
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
axis([0 800 time(1)/u.ps time(end)/u.ps])
caxis([-10*max(max(strainMap)) 10*max(max(strainMap))]);
colormap fireice(255);
colorbar('location', 'EastOutside');
xlabel('z [nm]');
ylabel('Delay [ps]');
title('Strain [%]');

%% Initilaize dynamic XRD Simulation
disp('Initilaize dynamic XRD Simulation');
D = XRDdyn(S,forceRecalc,E,sp); % set main parameters
D.setCacheDir(cacheDir); % set the cache directory
D.setQz(qz); % set q_z range

%% Homogeneous Dynamic XRD Simulations
% Here we calculate the rocking curve for no strain at all for the full 
% sample .
Rh = D.homogeneousReflectivity(); % thats all

% define an instrumental function to convolute the result with:
mu          = 0.1;
width       = 0.001*u.ang^-1;
instFunc    = @(qz)(pseudo_voigt(qz,width,mu));

% carry out the convolution with the instrumental function 
[Rhi xc] = D.convWithInstrumentFunction(Rh,qz,instFunc);

% Plot the Results of the homogeneous Dynamic XRD Simulations
figure(3)
semilogy(qz/u.ang^-1, Rh, '-b', 'LineWidth', 1);
hold on
semilogy(qz/u.ang^-1, instFunc(qz-xc)/2, '-k', 'LineWidth', 1);
semilogy(qz/u.ang^-1, Rhi, '-r', 'LineWidth', 2);
hold off
axis([qz(1)/u.ang^-1 qz(end)/u.ang^-1 1e-6 1]);
set(gca,'YScale', 'log');
grid on; box on;
title('Homogeneous XRD');
xlabel('q_z [Ang.^{-1}]'); 
ylabel('Reflectivity');
lh = legend('Rock. Curve', 'Inst. Func.', 'Conv. Rock. Curve');
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

% convolute with instrument function first
RI = D.convWithInstrumentFunction(R,qz,instFunc);

figure(4)
kk = surf(qz/u.ang^-1,time/u.ps,log10(RI));
set(kk, 'LineStyle', 'none');
axis([qz(1)/u.ang^-1 qz(end)/u.ang^-1 time(1)/u.ps time(end)/u.ps])
box on; colorbar; colormap(fireice(256));
caxis([-5 -2]);
xlabel('q_z [Ang^{-1}])'); 
ylabel('Delay [ps]');
title('Transient X-ray Reflectivity');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');

%%
tPlot   = (0:7.2:60)*u.ps;
figure(5)
kk = waterfall(qz/u.ang^-1,tPlot/u.ps,log10(RI(finderb(tPlot,time),:)));
set(kk, 'LineWidth', 2);
xlabel('q_z [Ang^{-1}])'); 
ylabel('Delay [ps]');
zlabel('Reflectivity');
title('Transient Substrate Side Peaks');
view([0,0]);
zlim([-4.5 -2.5]);
xlim([3.225 3.245]);
caxis([-4 -2]);
colormap(winter);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out');
box off; grid on;
%%%
disp('finished');

%% References
%
% # M. Herzog, A. Bojahr, J. Goldshteyn, W. Leitenberger, I. Vrejoiu, 
% D. Khakhulin, M. Wulff et al. (2012). _Detecting optically synthesized 
% quasi-monochromatic sub-terahertz phonon wavepackets by ultrafast x-ray 
% diffraction_. Applied Physics Letters, 100(9), 94101. 
% doi:10.1063/1.3688492