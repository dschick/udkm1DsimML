function y = gauss(x,varargin)
    % initialize input parser and define defaults and validators
    p = inputParser;
    p.addRequired('x'                   , @isnumeric);
    p.addOptional('s'           , 1     , @isnumeric);
    p.addOptional('widthType' , 'std' , @(x)(find(strcmp(x,{'std', 'FWHM', 'var', 'HWHM'}))));
    p.addOptional('normalize' , true  , @islogical);
    % parse the input
    p.parse(x,varargin{:});
    % assign parser results to object properties
    
    switch p.Results.widthType
        case 'FWHM' % Full Width at Half Maximum
            s = p.Results.s/(2*sqrt(2*log(2)));
        case 'var' % variance
            s = sqrt(p.Results.s);
        case 'HWHM' % Half Width at Half Maximum
            s = p.Results.s/(sqrt(2*log(2)));
        otherwise % standard derivation
            s = p.Results.s;
    end
    
    if p.Results.normalize == true
        a = 1/sqrt(2*pi*s^2); % normalize area to 1
    else
        a = 1; % normalize amplitude to 1
    end
        
    y = a * exp(-(x.^2)./(2*s^2));    
end