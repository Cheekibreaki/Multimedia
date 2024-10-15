function predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize, width, height)
    predictedFrame = zeros(height, width, 'uint8');
    [mvHeight, mvWidth, ~] = size(motionVectors);

    for y = 1:mvHeight
        for x = 1:mvWidth
            mv = motionVectors(y, x, :);
            refY = max(1, min(height-blockSize+1, (y-1)*blockSize+1 + mv(1)));
            refX = max(1, min(width-blockSize+1, (x-1)*blockSize+1 + mv(2)));
            predictedFrame((y-1)*blockSize+1:min(y*blockSize,height), ...
                           (x-1)*blockSize+1:min(x*blockSize,width)) = ...
                referenceFrame(refY:refY+blockSize-1, refX:refX+blockSize-1);
        end
    end
end