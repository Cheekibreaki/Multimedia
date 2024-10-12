function [motionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange)
    [height, width] = size(currentFrame);
    motionVectors = zeros(ceil(height/blockSize), ceil(width/blockSize), 2);
    totalMAE = 0;
    blockCount = 0;
    
    paddedReference = padarray(referenceFrame, [searchRange searchRange], 'replicate', 'both');
    
    % Pre-compute all possible blocks in the search range
    [dy, dx] = meshgrid(-searchRange:searchRange, -searchRange:searchRange);
    searchOffsets = [dy(:) dx(:)];
    
    for y = 1:blockSize:height
        for x = 1:blockSize:width
            blockCount = blockCount + 1;
            currentBlock = currentFrame(y:min(y+blockSize-1,height), x:min(x+blockSize-1,width));
            [actualBlockHeight, actualBlockWidth] = size(currentBlock);
            
            % Vectorized MAE calculation
            MAEs = zeros(size(searchOffsets, 1), 1);
            for i = 1:size(searchOffsets, 1)
                refY = y + searchRange + searchOffsets(i, 1);
                refX = x + searchRange + searchOffsets(i, 2);
                referenceBlock = paddedReference(refY:refY+actualBlockHeight-1, refX:refX+actualBlockWidth-1);
                MAEs(i) = sum(abs(double(currentBlock(:)) - double(referenceBlock(:)))) / (actualBlockHeight * actualBlockWidth);
            end
            
            % Find the best motion vector
            [bestMAE, bestIdx] = min(MAEs);
            bestMV = searchOffsets(bestIdx, :);
            
            motionVectors(ceil(y/blockSize), ceil(x/blockSize), :) = bestMV;
            totalMAE = totalMAE + bestMAE;
            
            % Early termination (optional)
            if bestMAE < 1  % You can adjust this threshold
                break;
            end
        end
    end
    
    avgMAE = totalMAE / blockCount;
end