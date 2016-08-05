function cellMatrix = exponentToCellMatrix(cellMatrix,exponent)
    % apply exponent to each matric in a cellArray of matrices
    if exponent > 1
        for i = 1:length(cellMatrix)
            cellMatrix{i} = cellMatrix{i}^exponent;
        end
    end
end

