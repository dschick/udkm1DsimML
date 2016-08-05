function R = multiplyCellMatrices(A, B)
    % multiply to 1D CellArrays of Matricies
    if length(A) ~= length(B)
        error('Both matrices must have the same length');
    end

    R = cell(length(A),1);
    for i = 1:length(A)
        R{i} = A{i} * B{i};
    end
end
