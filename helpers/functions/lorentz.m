function y = lorentz(x,w)
    % w is FWHM
    if (nargin == 1)
        w = 1;
    end
    
    y = 2/pi * w ./ (4*x.^2 + w^2);
    
end