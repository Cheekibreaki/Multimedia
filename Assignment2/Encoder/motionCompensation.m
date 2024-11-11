function predictedFrame = motionCompensation(referenceFrames, motionVectors, blockSize, width, height, FMEEnable)
    % Parameters:
    % referenceFrame - the reference frame (previous frame or hypothetical frame)
    % motionVectors  - motion vectors for each block
    % blockSize      - size of each block 
    
    
    % Initialize the predicted frame with zeros
    predictedFrame = zeros(height, width, 'double');
    
    
    % Loop over each block based on the pixel positions (row, col)
    for row = 1:blockSize:height
        for col = 1:blockSize:width
            % Convert row and col to block indices
            blockY = (row - 1) / blockSize + 1;  
            blockX = (col - 1) / blockSize + 1;  

            % Get the motion vector for this block (dy, dx)
            dy = motionVectors(blockY, blockX, 1);  % Vertical offset
            dx = motionVectors(blockY, blockX, 2);  % Horizontal offset
            refIdx = motionVectors(blockY, blockX, 3) + 1;  % Reference frame index

            % Extract the appropriate reference frame
            referenceFrame = referenceFrames{refIdx};
            
            % Compute the coordinates of the reference block based on the motion vector
            refRowStart = row + dy;
            refColStart = col + dx;

            % Extract the reference block from the reference frame
            if FMEEnable
                refRowStart = 2*row-1 + dy;
                refColStart = 2*col-1 + dx;
                refBlock = referenceFrame(refRowStart:2:(refRowStart + 2* blockSize - 2), refColStart:2:(refColStart + 2 * blockSize - 2));
            else
                refBlock = referenceFrame(refRowStart:(refRowStart + blockSize - 1), refColStart:(refColStart + blockSize - 1));
            end
            
            
            % Place the reference block into the predicted frame
            predictedFrame(row:(row + blockSize - 1), col:(col + blockSize - 1)) = double(max(0,min(255,refBlock)));
        end
    end
end
