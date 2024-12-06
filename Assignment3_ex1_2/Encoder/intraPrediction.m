function [approximatedPredictedFrame, predictionModes] = intraPrediction(currentFrame, blockSize,dct_blockSize,baseQP,RCflag,per_block_row_budget, bitCountPerRow)
    [height, width] = size(currentFrame);
    approximatedPredictedFrame = zeros(size(currentFrame), 'double');
    predictionModes = int32(zeros(ceil(height/blockSize), ceil(width/blockSize)));
    approximatedReconstructedFrame(1:height,1:width) = zeros(size(currentFrame), 'double');
    quantizedResiduals = zeros(size(currentFrame));
    row_bits_used = per_block_row_budget;
    for y = 1:blockSize:height
        if RCflag
                next_row_budget = per_block_row_budget + (per_block_row_budget - row_bits_used);
                baseQP = findCorrectQP(next_row_budget,bitCountPerRow);
                
            row_bits_used = 0;
        end
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
            
            
            quantizedBlock = quantization(residualBlock, dct_blockSize,blockSize,blockSize,baseQP);

            quantizedBlock1d = zigzag(quantizedBlock);
            residuesRLE = rle_encode(quantizedBlock1d); 
            encodedResidues = exp_golomb_encode(residuesRLE);
            
            if RCflag
                % Calculate bits used by this block
                row_bits_used = row_bits_used + length(encodedResidues);
            end

            approximatedresidualBlock = invquantization(quantizedBlock,dct_blockSize,blockSize,blockSize,baseQP);
            approximatedReconstructedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = double(max(0, min (255, approximatedresidualBlock + approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1))));      
        end
    end
end

function zigzag_order = zigzag(matrix)
    % Zigzag function to reorder a 2D matrix into a 1D array
    [rows, cols] = size(matrix);
    zigzag_order = zeros(1, rows * cols);
    index = 1;
    for s = 1:(rows + cols - 1)

        % Even diagonals: top-right to bottom-left
        for i = max(1, s - cols + 1):min(rows, s)
            j = s - i + 1;
            zigzag_order(index) = matrix(i, j);
            index = index + 1;
        end
    end
end




function encoded = exp_golomb_encode(data)
    % Encode data using Exponential-Golomb coding based on specific rules
    % and return a 1D array of binary values (0s and 1s)
    encoded = [];  % Initialize an empty array to store binary values
    
    for i = 1:length(data)
        x = data(i);
        % Map value to a non-negative integer
        if x > 0
            value = 2 * x;  % Positive x to odd integer (2x)
        else
            value = (-2 * x) + 1;  % Non-positive x to even integer (-2x + 1)
        end
        
        bin_code = dec2bin(value) - '0';  % Convert to binary and get an array of 1s and 0s
        leading_zeros = floor(log2(value)) + 1;  % Compute the number of leading zeros
        code = [zeros(1, leading_zeros - 1), 1, bin_code(2:end)];  % Assemble the code
        
        encoded = [encoded, code];  % Append binary code to the output array
    end
end




function encoded = rle_encode(data)
    % Encode data using a custom Run-Length Encoding (RLE) scheme
    % Return a 1D array containing the encoded values
    encoded = [];
    
    i = 1;
    while i <= length(data)
        if data(i) ~= 0
            % Count the number of consecutive non-zero values
            non_zero_count = 0;
            start_idx = i;
            while i <= length(data) && data(i) ~= 0
                non_zero_count = non_zero_count + 1;
                i = i + 1;
            end
            % Encode the non-zero run: prefix with -count followed by the values
            encoded = [encoded, -non_zero_count, data(start_idx:start_idx + non_zero_count - 1)];
        else
            % Count the number of consecutive zero values
            zero_count = 0;
            while i <= length(data) && data(i) == 0
                zero_count = zero_count + 1;
                i = i + 1;
            end

            if i == length(data)+1
                zero_count = 0;
            end
            % Encode the zero run: prefix with +count
            encoded = [encoded, zero_count];
        end
    end
    
    % Add a trailing 0 if the rest of the elements are zeros
end



