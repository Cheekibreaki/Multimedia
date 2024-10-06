function predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize)
    % Parameters:
    % referenceFrame - the reference frame (previous frame or hypothetical frame)
    % motionVectors  - motion vectors for each block
    % blockSize      - size of each block 
    
    % Get the size of the reference frame
    [height, width] = size(referenceFrame);
    
    % Initialize the predicted frame with zeros
    predictedFrame = zeros(height, width, 'uint8');
    
    
    % Loop over each block based on the pixel positions (row, col)
    for row = 1:blockSize:height
        for col = 1:blockSize:width
            % Convert row and col to block indices
            blockY = (row - 1) / blockSize + 1;  
            blockX = (col - 1) / blockSize + 1;  

            % Get the motion vector for this block (dy, dx)
            dy = motionVectors(blockY, blockX, 1);  % Vertical offset
            dx = motionVectors(blockY, blockX, 2);  % Horizontal offset
            
            % Compute the coordinates of the reference block based on the motion vector
            refRowStart = row + dy;
            refColStart = col + dx;

            % Extract the reference block from the reference frame
            refBlock = referenceFrame(refRowStart:(refRowStart + blockSize - 1), refColStart:(refColStart + blockSize - 1));
            
            % Place the reference block into the predicted frame
            predictedFrame(row:(row + blockSize - 1), col:(col + blockSize - 1)) = refBlock;
        end
    end
end
