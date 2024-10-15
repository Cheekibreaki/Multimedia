function [motionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange)
    [height, width] = size(currentFrame);
    motionVectors = zeros(ceil(height/blockSize), ceil(width/blockSize), 2);
    totalMAE = 0;
    blockCount = 0;

    paddedReference = padarray(referenceFrame, [searchRange searchRange], 'replicate', 'both');

    for y = 1:blockSize:height
        for x = 1:blockSize:width
            blockCount = blockCount + 1;
            currentBlock = currentFrame(y:min(y+blockSize-1,height), x:min(x+blockSize-1,width));
            [actualBlockHeight, actualBlockWidth] = size(currentBlock);
            
            bestMAE = inf;
            bestMV = [0, 0];
            
            for dy = -searchRange:searchRange
                for dx = -searchRange:searchRange
                    refY = y + searchRange + dy;
                    refX = x + searchRange + dx;
                    referenceBlock = paddedReference(refY:refY+actualBlockHeight-1, refX:refX+actualBlockWidth-1);
                    
                    MAE = sum(abs(double(currentBlock(:)) - double(referenceBlock(:)))) / (actualBlockHeight * actualBlockWidth);
                    
                    if MAE < bestMAE || (MAE == bestMAE && norm([dy, dx]) < norm(bestMV))
                        bestMAE = MAE;
                        bestMV = [dy, dx];
                    end
                end
            end
            
            motionVectors(ceil(y/blockSize), ceil(x/blockSize), :) = bestMV;
            totalMAE = totalMAE + bestMAE;
        end
    end

    avgMAE = totalMAE / blockCount;
end