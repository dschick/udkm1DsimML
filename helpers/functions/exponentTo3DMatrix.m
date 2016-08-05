function Matrix = exponentTo3DMatrix(Matrix,exponent)
    % apply exponent to each matric in a cellArray of matrices
    if exponent > 1
        for i = 1:size(Matrix,3)
            Matrix(:,:,i) = Matrix(:,:,i)^exponent;
        end
    end
end

