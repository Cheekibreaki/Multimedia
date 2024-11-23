function [approximatedPredictedFrame, predictionModes] = intraPrediction_block_parallel(currentFrame, blockSize, dct_blockSize, baseQP)
    % Type 2: Block-level parallelism (Wavefront)
    % Implements diagonal (wavefront) processing with parallelism using `parfor`.
    
    [height, width] = size(currentFrame);
    approximatedPredictedFrame = zeros(size(currentFrame), 'double');
    predictionModes = int32(zeros(ceil(height / blockSize), ceil(width / blockSize)));
    approximatedReconstructedFrame = zeros(size(currentFrame), 'double');
    
    % Number of blocks in each dimension
    numBlocksX = ceil(width / blockSize);
    numBlocksY = ceil(height / blockSize);

    % Loop over diagonals
    % Each wavefront can be encoded in parallel
    for diagonal = 1:(numBlocksX + numBlocksY - 1)
        % Determine the range of blocks in this diagonal
        startBlockY = max(1, diagonal - numBlocksX + 1);
        endBlockY = min(diagonal, numBlocksY);

        % temp storage for results in this particular wavefront 
        tempPredictedFrame = zeros(size(currentFrame), 'double');
        tempReconstructedFrame = zeros(size(currentFrame), 'double');
        tempPredictionModes = zeros(numBlocksY, numBlocksX, 'int32');
        
        % Process the current diagonal using `parfor`
        parfor blockIdx = startBlockY:endBlockY
            % Create independent local variables for each `parfor` iteration
            localPredictedFrame = zeros(size(currentFrame), 'double');
            localReconstructedFrame = zeros(size(currentFrame), 'double');
            localPredictionModes = zeros(numBlocksY, numBlocksX, 'int32');

            blockY = blockIdx;
            blockX = diagonal - blockY + 1;

            % Calculate the pixel range for the current block
            x = (blockX - 1) * blockSize + 1;
            y = (blockY - 1) * blockSize + 1;
            actualBlockWidth = min(blockSize, width - x + 1);
            actualBlockHeight = min(blockSize, height - y + 1);

            % Block processing
            if y == 1 && x == 1
                % For the first block, predict with mid-gray
                localPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = 128;
                localPredictionModes(blockY, blockX) = 2; % 2 for mid-gray prediction
            else
                % Horizontal prediction
                if x > 1
                    horizPred = repmat(approximatedReconstructedFrame(y:y+actualBlockHeight-1, x-1), 1, actualBlockWidth);
                else
                    horizPred = repmat(128, actualBlockHeight, actualBlockWidth);
                end

                % Vertical prediction
                if y > 1
                    vertPred = repmat(approximatedReconstructedFrame(y-1, x:x+actualBlockWidth-1), actualBlockHeight, 1);
                else
                    vertPred = repmat(128, actualBlockHeight, actualBlockWidth);
                end

                % Calculate MAE for both predictions
                currentBlock = currentFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1);
                maeHoriz = mean(abs(double(currentBlock) - double(horizPred)), 'all');
                maeVert = mean(abs(double(currentBlock) - double(vertPred)), 'all');
                if maeHoriz <= maeVert
                    localPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = horizPred;
                    localPredictionModes(blockY, blockX) = 0; % 0 for horizontal
                else
                    localPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = vertPred;
                    localPredictionModes(blockY, blockX) = 1; % 1 for vertical
                end
            end

            % Compute residual, quantize, reconstruct
            residualBlock = currentFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) - localPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1);
            quantizedBlock = quantization(residualBlock, dct_blockSize, blockSize, blockSize, baseQP);
            approximatedResidualBlock = invquantization(quantizedBlock, dct_blockSize, blockSize, blockSize, baseQP);
            localReconstructedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = double(max(0, min(255, approximatedResidualBlock + localPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1))));

            % Aggregate local results into the temporary arrays for the wavefront
            tempPredictedFrame = tempPredictedFrame + localPredictedFrame;
            tempReconstructedFrame = tempReconstructedFrame + localReconstructedFrame;
            tempPredictionModes = tempPredictionModes + localPredictionModes;
        end

        % Update shared variables after processing the diagonal
        approximatedPredictedFrame = approximatedPredictedFrame + tempPredictedFrame;
        approximatedReconstructedFrame = approximatedReconstructedFrame + tempReconstructedFrame;
        predictionModes = predictionModes + tempPredictionModes;
    end
end
