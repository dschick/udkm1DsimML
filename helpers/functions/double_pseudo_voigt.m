function y = double_pseudo_voigt(x,w,mu,a,c)
    if (nargin < 2)
        w = 1;
    end
    if (nargin < 3)
        mu = 0.5;
    end
    if (nargin < 4)
        a = 1;
    end
    if (nargin < 5)
        c = 1;
    end
    
    y = pseudo_voigt(x,w,mu) + a *pseudo_voigt(x-c,w,mu);
    
end