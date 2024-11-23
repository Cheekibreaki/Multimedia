function originalValues = diffDecoding(values, type)
    % Generalized differential decoding function for both prediction modes and motion vectors
    % Input:
    %   diffValues - Matrix of differential values (2D for modes, 3D for motion vectors)
    %   type       - 'modes' for prediction modes, 'mv' for motion vectors
    % Output:
    %   originalValues - Decoded matrix with the original values

    [numBlocksY, numBlocksX, numDims] = size(values);  % Get dimensions
    originalValues = zeros(size(values));  % Initialize output matrix

    for i = 1:numBlocksY
        for j = 1:numBlocksX
            if j == 1  % First block in the row
                if strcmp(type, 'modes')
                    % Assume previous mode is Horizontal (0)
                    originalValues(i, j) = values(i, j) + 0;
                elseif strcmp(type, 'mv')
                    % Assume previous MV is (0, 0, 0)
                    originalValues(i, j, 1) = values(i, j, 1) + 0;
                    originalValues(i, j, 2) = values(i, j, 2) + 0;
                    originalValues(i, j, 3) = values(i, j, 3) + 0;
                end
            else  % Add the differential value to the previous block's value
                if strcmp(type, 'modes')
                    originalValues(i, j) = values(i, j) + originalValues(i, j - 1);
                elseif strcmp(type, 'mv')
                    originalValues(i, j, 1) = values(i, j, 1) + originalValues(i, j - 1, 1);
                    originalValues(i, j, 2) = values(i, j, 2) + originalValues(i, j - 1, 2);
                    originalValues(i, j, 3) = values(i, j, 3) + originalValues(i, j - 1, 3);
                end
            end
        end
    end
end
