function predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize)
    [height, width] = size(referenceFrame);
    [numBlocksY, numBlocksX, ~] = size(motionVectors);
    
    % Pre-compute block indices
    [blockIndicesY, blockIndicesX] = ndgrid(1:blockSize, 1:blockSize);
    
    predictedFrame = zeros(size(referenceFrame), 'uint8');
    
    for blockY = 1:numBlocksY
        for blockX = 1:numBlocksX
            y = (blockY - 1) * blockSize + 1;
            x = (blockX - 1) * blockSize + 1;
            
            mvY = motionVectors(blockY, blockX, 1);
            mvX = motionVectors(blockY, blockX, 2);
            
            refY = max(1, min(y + mvY, height - blockSize + 1));
            refX = max(1, min(x + mvX, width - blockSize + 1));
            
            predictedBlock = referenceFrame(refY:min(refY+blockSize-1,height), refX:min(refX+blockSize-1,width));
            predictedFrame(y:min(y+blockSize-1,height), x:min(x+blockSize-1,width)) = predictedBlock;
        end
    end
end