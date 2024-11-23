function [diffValues,last_value_output] = diffEncoding_block(curr_value,type,last_value)
    % Generalized differential encoding function for both prediction modes and motion vectors
    % Input:
    %   curr_value - 2D (for prediction modes) or 3D (for motion vectors) matrix
    %   type   - 'modes' for prediction modes, 'mv' for motion vectors
    % Output:
    %   diffValues - Matrix with differential curr_value of the same size as inputs

    [numBlocksY, numBlocksX, numDims] = size(curr_value);  % Get dimensions
    diffValues = zeros(size(curr_value));  % Initialize output matrix

    for i = 1:numBlocksY
        for j = 1:numBlocksX
            if j == 1 && i == 1 % First block in the row
                if strcmp(type, 'modes')
                    % Assume previous mode is Horizontal (0)
                    diffValues(i, j) = curr_value(i, j) - last_value(i,j);  
                elseif strcmp(type, 'mv')
                    % Assume previous MV is (0, 0, 0)
                    diffValues(i, j, 1) = curr_value(i, j, 1) - last_value(i,j,1);
                    diffValues(i, j, 2) = curr_value(i, j, 2) - last_value(i,j,2);
                    diffValues(i, j, 3) = curr_value(i, j, 3) - last_value(i,j,3); % Reference frame index
                end
            elseif j == 1 % First block in the column, use previous row's last value
                if strcmp(type, 'modes')
                    diffValues(i, j) = curr_value(i, j) - curr_value(i-1, numBlocksX);
                elseif strcmp(type, 'mv')
                    diffValues(i, j, 1) = curr_value(i, j, 1) - curr_value(i-1, numBlocksX, 1);
                    diffValues(i, j, 2) = curr_value(i, j, 2) - curr_value(i-1, numBlocksX, 2);
                    diffValues(i, j, 3) = curr_value(i, j, 3) - curr_value(i-1, numBlocksX, 3); % Reference frame index
                end
            else % Compare with the left neighbor
                if strcmp(type, 'modes')
                    diffValues(i, j) = curr_value(i, j) - curr_value(i, j - 1);
                elseif strcmp(type, 'mv')
                    diffValues(i, j, 1) = curr_value(i, j, 1) - curr_value(i, j - 1, 1);
                    diffValues(i, j, 2) = curr_value(i, j, 2) - curr_value(i, j - 1, 2);
                    diffValues(i, j, 3) = curr_value(i, j, 3) - curr_value(i, j - 1, 3);  % Reference frame index
                end
            end
        end
    end
    last_value_output = curr_value(numBlocksY,numBlocksX,:);
end