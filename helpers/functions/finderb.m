%% finderb
% Binary search algorithm for sorted lists.
% Searches for the first index _i_ of list where _key_ <= _list(i)_.
% _key_ can be a scalar or a vector of keys. 
% _list_ must be a sorted vector.
% author: André Bojahr
% licence: BSD
function i = finderb(key,list)
    n = length(key);
    i = zeros(1,n);

    if n > 500000; % if t is too long, we parallize it
        parfor (m = 1:n ,4)
            i(m) = finderb2(key(m),list);
        end%parfor
    else
        for m = 1:n
            i(m) = finderb2(key(m),list);
        end%for
    end%if
end%function

%% nested subfunction
function i = finderb2(key,list)
    a = 1;              % start of intervall
    b = length(list);   % end of intervall    
    
    % if the key is smaller than the first element of the
    % list we return 1
    if key < list(1)
        i = 1;
        return;
    end%if
    
    while (b-a) > 1 % loop until the intervall is larger than 1
        c = floor((a+b)/2); % center of intervall
        if key < list(c)
            % the key is in the left half-intervall
            b = c;
        else
            % the key is in the right half-intervall
            a = c;
        end%if        
    end%while
    
    i = a;
end%function