function diffValues = diffEncoding(values,type)
    % Generalized differential encoding function for both prediction modes and motion vectors
    % Input:
    %   values - 2D (for prediction modes) or 3D (for motion vectors) matrix
    %   type   - 'modes' for prediction modes, 'mv' for motion vectors
    % Output:
    %   diffValues - Matrix with differential values of the same size as input

    [numBlocksY, numBlocksX, numDims] = size(values);  % Get dimensions
    diffValues = zeros(size(values));  % Initialize output matrix

    for i = 1:numBlocksY
        for j = 1:numBlocksX
            if j == 1  % First block in the row
                if strcmp(type, 'modes')
                    % Assume previous mode is Horizontal (0)
                    diffValues(i, j) = values(i, j) - 0;  
                elseif strcmp(type, 'mv')
                    % Assume previous MV is (0, 0, 0)
                    diffValues(i, j, 1) = values(i, j, 1) - 0; 
                    diffValues(i, j, 2) = values(i, j, 2) - 0;
                    diffValues(i, j, 3) = values(i, j, 3) - 0;  % Reference frame index
                end
            else  % Compare with the left neighbor
                if strcmp(type, 'modes')
                    diffValues(i, j) = values(i, j) - values(i, j - 1);
                elseif strcmp(type, 'mv')
                    diffValues(i, j, 1) = values(i, j, 1) - values(i, j - 1, 1);
                    diffValues(i, j, 2) = values(i, j, 2) - values(i, j - 1, 2);
                    diffValues(i, j, 3) = values(i, j, 3) - values(i, j - 1, 3);  % Reference frame index
                end
            end
        end
    end
end