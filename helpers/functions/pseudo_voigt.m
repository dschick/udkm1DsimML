function y = pseudo_voigt(x,w,mu)
    if (nargin == 1)
        w = 1;
        mu = 0.5;
    elseif (nargin == 2)
        mu = 0.5; 
    end
    
    y = mu * 1 ./ (1 + (x./w).^2) + (1-mu) .* exp(-log(2) .* (x./w).^2);
    
end