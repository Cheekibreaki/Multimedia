% function [predictedFrame, predictionModes] = intraPrediction(currentFrame, blockSize)
%     [height, width] = size(currentFrame);
%     predictedFrame = zeros(size(currentFrame), 'uint8');
%     predictionModes = int32(zeros(ceil(height/blockSize), ceil(width/blockSize)));
%     predictionBlock = zeros(size(blockSize),'double');
% 
% 
%     for y = 1:blockSize:height
%         for x = 1:blockSize:width
%             actualBlockHeight = min(blockSize, height-y+1);
%             actualBlockWidth = min(blockSize, width-x+1);
% 
%             if y == 1 && x == 1
%                 % For the first block, predict with mid-gray
%                 predictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = 128;
%                 predictionModes(1, 1) = 2; % 2 for mid-gray prediction
%             else
%                 % Horizontal prediction
%                 if x > 1
%                     horizPred = repmat(currentFrame(y:y+actualBlockHeight-1, x-1), 1, actualBlockWidth);
%                 else
%                     horizPred = repmat(128, actualBlockHeight, actualBlockWidth);
%                 end
% 
%                 % Vertical prediction
%                 if y > 1
%                     vertPred = repmat(currentFrame(y-1, x:x+actualBlockWidth-1), actualBlockHeight, 1);
%                 else
%                     vertPred = repmat(128, actualBlockHeight, actualBlockWidth);
%                 end
% 
%                 % Calculate MAE for both predictions
%                 currentBlock = currentFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1);
%                 maeHoriz = mean(abs(double(currentBlock) - double(horizPred)), 'all');
%                 maeVert = mean(abs(double(currentBlock) - double(vertPred)), 'all');
% 
%                 % Choose the best prediction
%                 if maeHoriz <= maeVert
%                     predictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = horizPred;
%                     predictionModes(ceil(y/blockSize), ceil(x/blockSize)) = 0; % 0 for horizontal
%                 else
%                     predictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = vertPred;
%                     predictionModes(ceil(y/blockSize), ceil(x/blockSize)) = 1; % 1 for vertical
%                 end
%             end
%         end
%     end
% end



function [approximatedPredictedFrame, predictionModes,approximatedReconstructedFrame] = intraPrediction(currentFrame, blockSize,dct_blockSize,baseQP)
    [height, width] = size(currentFrame);
    approximatedPredictedFrame = zeros(size(currentFrame), 'double');
    predictionModes = int32(zeros(ceil(height/blockSize), ceil(width/blockSize)));


    predictionBlock = zeros(size(blockSize),'double');
    approximatedReconstructedFrame(1:height,1:width) = zeros(size(currentFrame), 'double');
    
    
    for y = 1:blockSize:height
        for x = 1:blockSize:width
            actualBlockHeight = min(blockSize, height-y+1);
            actualBlockWidth = min(blockSize, width-x+1);
            
            if y == 1 && x == 1
                % For the first block, predict with mid-gray
                approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = 128;
                predictionModes(1, 1) = 2; % 2 for mid-gray prediction
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
                    approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = horizPred;
                    predictionModes(ceil(y/blockSize), ceil(x/blockSize)) = 0; % 0 for horizontal
                else
                    approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = vertPred;
                    predictionModes(ceil(y/blockSize), ceil(x/blockSize)) = 1; % 1 for vertical
                end
            end
            
            residualBlock = currentFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) - approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1);
            quantizedBlock = quantization(residualBlock, dct_blockSize,blockSize,blockSize,baseQP)
            approximatedresidualBlock = invquantization(quantizedBlock,dct_blockSize,blockSize,blockSize,baseQP)
            approximatedReconstructedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = approximatedresidualBlock + approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1);
        end
    end
end