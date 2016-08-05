classdef constants
    properties(Constant)
        c_0     = 2.99792458e8;     % [m/s]
        e_0     = 1.60217733e-19;   % [C] 
        h       = 6.62606896e-34;   % [Js[
        hbar    = 1.05457266e-34;   % [Js[
        eps_0   = 8.854187871e-12;  % [A^2 s^4 kg^-1 m^-3]
        m_0     = 9.1093897e-31;    % [kg]
        r_e     = constants.e_0^2 / (4 * pi * constants.eps_0 * constants.m_0 * constants.c_0^2); % [m]
        u       = 1.660538921e-27;  % [kg]
        Na      = 6.02e23;          % [mol^-1]
    end
end