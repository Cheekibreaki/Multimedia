function [reconstructed_values,last_value] = diffDecoding_block(diffValues, type, last_value)
    % Generalized differential decoding function for both prediction modes and motion vectors
    % Input:
    %   diffValues - Matrix with differential values of the same size as the original inputs
    %   type       - 'modes' for prediction modes, 'mv' for motion vectors
    %   last_value - 2D (for prediction modes) or 3D (for motion vectors) matrix representing
    %                the previous value used for the first block
    % Output:
    %   reconstructed_values - Matrix with reconstructed values of the same size as inputs

    [numBlocksY, numBlocksX, numDims] = size(diffValues);  % Get dimensions
    reconstructed_values = zeros(size(diffValues));  % Initialize output matrix

    for i = 1:numBlocksY
        for j = 1:numBlocksX
            if j == 1 && i == 1 % First block in the row
                if strcmp(type, 'modes')
                    % Assume previous mode is Horizontal (0)
                    reconstructed_values(i, j) = diffValues(i, j) + last_value(i, j);
                elseif strcmp(type, 'mv')
                    % Assume previous MV is (0, 0, 0)
                    reconstructed_values(i, j, 1) = diffValues(i, j, 1) + last_value(i, j, 1);
                    reconstructed_values(i, j, 2) = diffValues(i, j, 2) + last_value(i, j, 2);
                    reconstructed_values(i, j, 3) = diffValues(i, j, 3) + last_value(i, j, 3); % Reference frame index
                end
            elseif j == 1 % First block in the column, use previous row's last value
                if strcmp(type, 'modes')
                    reconstructed_values(i, j) = diffValues(i, j) + reconstructed_values(i - 1, numBlocksX);
                elseif strcmp(type, 'mv')
                    reconstructed_values(i, j, 1) = diffValues(i, j, 1) + reconstructed_values(i - 1, numBlocksX, 1);
                    reconstructed_values(i, j, 2) = diffValues(i, j, 2) + reconstructed_values(i - 1, numBlocksX, 2);
                    reconstructed_values(i, j, 3) = diffValues(i, j, 3) + reconstructed_values(i - 1, numBlocksX, 3); % Reference frame index
                end
            else % Compare with the left neighbor
                if strcmp(type, 'modes')
                    reconstructed_values(i, j) = diffValues(i, j) + reconstructed_values(i, j - 1);
                elseif strcmp(type, 'mv')
                    reconstructed_values(i, j, 1) = diffValues(i, j, 1) + reconstructed_values(i, j - 1, 1);
                    reconstructed_values(i, j, 2) = diffValues(i, j, 2) + reconstructed_values(i, j - 1, 2);
                    reconstructed_values(i, j, 3) = diffValues(i, j, 3) + reconstructed_values(i, j - 1, 3); % Reference frame index
                end
            end
        end
    end

    last_value = reconstructed_values(numBlocksY,numBlocksX,:);
