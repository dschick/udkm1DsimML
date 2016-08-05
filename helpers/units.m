classdef units
% The units.m function returns a struct containing the SI values of 
    % many common units. The example below demonstrates how to use 
    % this struct to cleanly and effeciently perform matlab calculations 
    % using commonly encountered units and physical constants. The code 
    % is easily modified to include any non-standard unit or constant desired.
    %
    % To determine the exact syntax of a particular unit you would like to use,
    % simply run the function with no semicolon and all the units will be
    % displayed.
    %
    %  Using the units struct:
    %  --------------------------------------------------------
    %    First create a units struct by including in your code the following
    %    line
    %          u = units;
    %    Then:
    %    
    %          %To enter a number in a given unit, MULTIPLY by the unit:
    %                L = 5*u.in   % matlab automatically displays L in SI units
    %  
    %          % To display in a desired unit, simply divide by that unit
    %                L/u.ft       % displays L in ft.
    %  
    %          % To convert between units, MULTIPLY by the starting unit, and
    %          % DIVIDE by the ending unit:
    %                u.mi^2/u.ft^2  %displays the number of feet^2 in one mile^2
    %  
    %          %More complicated units can be obtained through arithmatic
    %                mach1 = 340.29*u.m/u.s;  %speed of sound 
    %  
    %             %Note... to make the speed of sound available wherever your
    %             %units struct is defined, simply write:
    %                u.mach1 = 340.29*u.m/u.s;   %mach1 now part of units struct
    % 
    %
    %
    % %------  BEGIN EXAMPLE CODE --------------------------------
    % %This is an example calculation that uses the units mfile to calculate the
    % %pressure at the bottom of a long vertically oriented pipe that is capped 
    % %at the bottom and filled with oil.
    % 
    % u = units;
    % pipeInnerDiameter = 4*u.in;     %4 inch inner diameter
    % pipeHeight = 30*u.ft;           %pipe sticks 30 feet up into the air
    % densityOfOil = 0.926*u.gm/u.cc; %density of oil as found on some random web site = .926 gm/cc
    % pipeCrossSectionArea = pi*(pipeInnerDiameter/2)^2;  %cross sectional area of pipe bore
    % volumeOfOil = pipeCrossSectionArea * pipeHeight;    %volume of oil that the pipe can hold
    % pressurePipeBottom = densityOfOil * u.g * pipeHeight;  %pressure formula from physics: P = rho*g*h.
    % forceOnPipeBottom = pressurePipeBottom * pipeCrossSectionArea;  %force exterted on bottom cap of the pipe.
    % 
    % %Note that each variable holds its value as expressed in SI units.  To
    % %express these values in different units, simply divide by the desired
    % %unit as shown below.
    % line1 = sprintf('A %2.3g inch diameter pipe sticking %3.3g meters into the air and filled',pipeInnerDiameter/u.in, pipeHeight/u.m);
    % line2 = sprintf('with %3.3g fluid ounces of oil will have a pressure at the bottom of %4.4g psi.',volumeOfOil/u.floz, pressurePipeBottom/u.psi);
    % line3 = sprintf('This will cause a total force of %5.5g lbs to press on the bottom cap of the pipe.',forceOnPipeBottom/u.lbf);
    % 
    % textVal = sprintf('\n\n%s\n%s\n%s\n',line1,line2,line3);
    % disp(textVal);
    %%------  END EXAMPLE CODE --------------------------------

    properties(Constant)

        %============ START THE ACTUAL CODE TO DEFINE THE UNITS STRUCT =========
        %-------- UNITS ------------------------------
        %------- length ----
        m = 1;
        km = 1e3*units.m;
        cm = 1e-2*units.m;
        mm = 1e-3*units.m;
        um = 1e-6*units.m;
        nm = 1e-9*units.m;
        ang = 1e-10*units.m;
        in = 2.54*units.cm;
        mil = 1e-3*units.in;
        ft = 12*units.in;
        yd = 3*units.ft;
        mi = 5280*units.ft;
        a0 = .529e-10*units.m;

        %------- Volume -------
        cc = (units.cm)^3;
        L = 1000*units.cc;
        mL = units.cc;
        floz = 29.5735297*units.cc;
        pint = 473.176475*units.cc;
        quart = 946.35295*units.cc;
        gal = 3.78541197*units.L;

        %----- mass ---------
        kg = 1;
        gm = 1e-3*units.kg;
        mg = 1e-3*units.gm;
        lb = 0.45359237*units.kg;
        oz = (1/16)*units.lb;
        amu = 1.66e-27*units.kg;

        %---- time -------
        s = 1;
        ms = 1e-3*units.s;
        us = 1e-6*units.s;
        ns = 1e-9*units.s;
        ps = 1e-12*units.s;
		fs = 1e-15*units.s;
		as = 1e-18*units.s;
        min = 60*units.s;
        hr = 60*units.min;
        day = 24*units.hr;
        yr = 365.242199*units.day; 

        %---- frequency ---- 
        Hz = 1/units.s;
        kHz = 1e3 *units.Hz;
        MHz = 1e6 *units.Hz;
        GHz = 1e9 *units.Hz;
        THz = 1e12 *units.Hz;

        %---- force -------
        N = 1;
        dyne = 1e-5*units.N;
        lbf = 4.44822*units.N;


        %----- energy -----
        J = 1;
        MJ = 1e6*units.J;
        kJ = 1e3*units.J;
        mJ = 1e-3*units.J;
        uJ = 1e-6*units.J;
        nJ = 1e-9*units.J;
        eV = 1.6022e-19*units.J;
        BTU = 1.0550559e3*units.J;
        kWh = 3.6e6*units.J;
        cal = 4.1868*units.J;
        kCal = 1e3*units.cal;

        %---- temperature ---
        K = 1;
        mK = 1e-3*units.K;
        uK = 1e-6*units.K;
        nK = 1e-9*units.K;

        %---- pressure -----
        Pa = 1;
        torr = 133.322*units.Pa;
        mtorr = 1e-3*units.torr;
        bar = 1e5*units.Pa;
        mbar = 1e-3*units.bar;
        atm = 1.013e5*units.Pa;
        psi = 6.895e3*units.Pa;



        %----- power --- ---
        W = 1;
        MW = 1e6*units.W;
        kW = 1e3*units.W;
        mW = 1e-3*units.W;
        uW = 1e-6*units.W;
        nW = 1e-9*units.W;
        pW = 1e-12*units.W;
        hp = 745.69987*units.W;

        %------ charge ------
        coul = 1;
        e = 1.6022e-19*units.coul;


        %------ Voltage -----
        V = 1;
        kV = 1e3*units.V;
        mV = 1e-3*units.V;
        uV = 1e-6*units.V;

        %----- Current ------
        A = 1;
        mA = 1e-3*units.A;
        uA = 1e-6*units.A;
        nA = 1e-9*units.A;

        %----magnetic field -----
        T = 1;
        gauss = 1e-4*units.T;
        
        % angles
        rad = 1;
        mrad = 1e-3*units.rad;
        deg = pi/180*units.rad;


        % constants are defined in constants class already in SI units
        % %----fundamental constants ----
        % g = 9.80665*units.m/units.s^2;
        % kB = 1.38e-23*units.J/units.K;
        % sigma_SB = 5.670e-8 * units.W/(units.m^2 * units.K^4);
        % h = 6.626e-34 * units.J*units.s;
        % hbar = units.h/(2*pi);
        % mu_B = 9.274e-24 * units.J/units.T;
        % mu_N = 5.0507866e-27 * units.J/units.T;
        % c = 2.99792458e8*units.m/units.s;
        % eps0 = 8.8541878176204e-12* units.coul/(units.V*units.m);
        % mu0 = 1.2566370614359e-6 * units.J/(units.m*units.A^2);
    end
end





