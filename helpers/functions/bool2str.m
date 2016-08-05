%% bool2str
% Returns the according string for a boolean input.
function str = bool2str(bool)
    if bool
        str = 'true';
    else
        str = 'false';
    end%if
end%function

