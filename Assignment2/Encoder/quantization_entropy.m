function [encodedMotionVector,encodedPredicitonModes,encodedResidues,quantizedResidues] = quantization_entropy(frame_type, motionVector3d, predicitonModes2d,residuals, RCflag,per_block_row_budget, bitCountPerRow,dct_blockSize,width,height,baseQP, vbs_matrix)
    % Input:
    % frame_type: 1 for I-frame, 0 for P-frame 
    % pred_diff: array containing differential prediction information (modes for intra or motion vectors for inter)
    % dct_coeffs: block of quantized DCT coefficients
    % block_size: size of the block (e.g., 4 for 4x4 blocks)

    % Frame Type Marker. Create the header to store frame type 
    frameTypeHeader = frame_type;  % 1 for I-frame, 0 for P-frame
    encodedMotionVector = [];
    encodedPredicitonModes = [];
    encodedResidues = [];
    if frame_type == 0  % Assuming data handling for P-frames
       if exist('vbs_matrix', 'var')
            % motionVector3d is 36x44x3, vbs_matrix is 36x44
            % Fetch a 2x2 block from each within vbs_matrix and motionVector3d
            
            [rows, cols] = size(vbs_matrix);
            resultmotionVector1d = [];
            
            % Loop through vbs_matrix and motionVector3d with a 2x2 block size

            previous_motion_vector_block = zeros(1, 1, 3);

            for row = 1:2:rows
                for col = 1:2:cols
                    % Extract the current 2x2 block from vbs_matrix
                    vbs_block = vbs_matrix(row:row+1, col:col+1);
                    
                    % Extract the corresponding 2x2x3 block from motionVector3d
                    motion_block = motionVector3d(row:row+1, col:col+1, :);
                    [motion_block,previous_motion_vector_block] = diffEncoding_block(motion_block,'mv',previous_motion_vector_block);
                    
                    if all(vbs_block(:) == 1)  % If this 2x2 block in vbs_matrix is all zeros
                        % Fetch the entire 2x2x3 block from motionVector3d
                        motionVector2d = reshape_3d_to_2d(motion_block);
                        motionVector1d = zigzag(motionVector2d);
                        motionVector1d = [1, motionVector1d];  % Prefix with 0 to indicate all zeros
                        resultmotionVector1d = [resultmotionVector1d, motionVector1d];  % Append to the result
                    elseif all(vbs_block(:) == 0)  % If this 2x2 block in vbs_matrix is all ones
                        % Fetch the 1x1x3 top-left element of the 2x2x3 motionVector3d block
                        top_left_block = motion_block(1, 1, :);
                        motionVector2d = reshape_3d_to_2d(top_left_block);
                        motionVector1d = zigzag(motionVector2d);
                        motionVector1d = [0, motionVector1d];  % Prefix with 1 to indicate all ones
                        resultmotionVector1d = [resultmotionVector1d, motionVector1d];  % Append to the result
                    else
                        error('ErrorID:1', 'a block can not contain both 1 and 0!');
                        %This shouldn't happen! something goes wrong
                    end
                end
            end
    
            % Encode the result using Exp-Golomb encoding
            encodedMotionVector = exp_golomb_encode(resultmotionVector1d);
       else


            % For block encoding if vbs_matrix does not exist
            motionVector2d = reshape_3d_to_2d(motionVector3d);
            motionVector1d = zigzag(motionVector2d);
            encodedMotionVector = exp_golomb_encode(motionVector1d);
        end

    elseif frame_type == 1

        if exist('vbs_matrix', 'var')
            [rows, cols] = size(vbs_matrix);
            resultpredictionMode_block_1d = [];
            
            % Loop through vbs_matrix and motionVector3d with a 2x2 block size

            previous_prediction_mode_block = 0;

            
            for row = 1:2:rows
                for col = 1:2:cols
                    % Extract the current 2x2 block from vbs_matrix
                    vbs_block = vbs_matrix(row:row+1, col:col+1);
                    
                    % Extract the corresponding 2x2x3 block from motionVector3d
                    predictionMode_block = predicitonModes2d(row:row+1, col:col+1, :);
                    [predictionMode_block,previous_prediction_mode_block] = diffEncoding_block(predictionMode_block,'modes',previous_prediction_mode_block);
                    
                    if all(vbs_block(:) == 1)  % If this 2x2 block in vbs_matrix is all ones
                        % Fetch the entire 2x2x3 block from motionVector3d
                        predictionMode_block_1d = zigzag(predictionMode_block);
                        predictionMode_block_1d = [1, predictionMode_block_1d];  % Prefix with 0 to indicate all zeros
                        resultpredictionMode_block_1d = [resultpredictionMode_block_1d, predictionMode_block_1d];  % Append to the result
                    elseif all(vbs_block(:) == 0)  % If this 2x2 block in vbs_matrix is all zeros
                        % Fetch the 1x1x3 top-left element of the 2x2x3 motionVector3d block
                        top_left_block = predictionMode_block(1, 1, :);
                        predictionMode_block_1d = zigzag(top_left_block);
                        predictionMode_block_1d = [0, predictionMode_block_1d];  % Prefix with 1 to indicate all ones
                        resultpredictionMode_block_1d = [resultpredictionMode_block_1d, predictionMode_block_1d];  % Append to the result
                    else
                        error('ErrorID:1', 'a block can not contain both 1 and 0!');
                        %This shouldn't happen! something goes wrong
                    end
                end
            end
    
            % Encode the result using Exp-Golomb encoding
            encodedPredicitonModes = exp_golomb_encode(resultpredictionMode_block_1d);
       else
            predicitonModes1d = zigzag(predicitonModes2d);
            encodedPredicitonModes = exp_golomb_encode(predicitonModes1d); 
        end
    end
    if(RCflag == false) 
        residues2d = zeros(size(residuals));
        for row = 1:dct_blockSize:height
            for col = 1:dct_blockSize:width
                block = residuals(row:row+dct_blockSize-1, col:col+dct_blockSize-1);
                dctBlock = dct2(double(block));
        
                % Quantization
                % QP = adaptiveQP(dct_blockSize, baseQP);
                Q = createQMatrix(size(block), baseQP);
                quantizedBlock = round(dctBlock ./ Q);
                residues2d(row:row+dct_blockSize-1, col:col+dct_blockSize-1) = quantizedBlock;
            end
        end
        quantizedResidues = residues2d;
        residues1d = zigzag(residues2d);
        residuesRLE = rle_encode(residues1d); 
        encodedResidues = exp_golomb_encode(residuesRLE);
    else
        final_encodedResidues = [];
        % Define maximum QP to prevent infinite loops
        maxQP = 10;  % Maximum QP value (standard range is 0 to 51)
        % Initialize a matrix to store quantized residues for reconstruction
        quantizedResidues = zeros(size(residuals));
        % Loop over the image in block rows
        for row = 1:dct_blockSize:height
            curr_per_block_row_budget = per_block_row_budget;
            % Loop over the blocks in the current block row
            for col = 1:dct_blockSize:width
                % Extract the current block
                block_row_end = min(row + dct_blockSize - 1, height);
                block_col_end = min(col + dct_blockSize -1, width);
                block = residuals(row:block_row_end, col:block_col_end);
                dctBlock = dct2(double(block));
    
                % Initialize QP with baseQP
                QP = baseQP;
                % Initialize variables
                encodedResidues = [];
                bits_used = 0;
                % Rate Control Loop: Adjust QP to meet the per-block row budget
                
                % Quantization with current QP
                Q = createQMatrix(size(block), QP);
                quantizedBlock = round(dctBlock ./ Q);

                % Encode the quantized block
                % Flatten the quantized block using zigzag scan
                quantizedBlock1d = zigzag(quantizedBlock);
                residuesRLE = rle_encode(quantizedBlock1d); 
                encodedResidues = exp_golomb_encode(residuesRLE);

                % Calculate bits used by this block
                bits_used = length(encodedResidues);

                % Check if bits used exceeds the remaining budget
                if bits_used <= curr_per_block_row_budget
                    % Sufficient budget; proceed with this QP
                    break;
                else
                    % Insufficient budget; increase QP to reduce bits
                        warning('Reached maximum QP at position (%d, %d); meet bit budget.', row, col);
                        break;
                    end
                end
                
                % Update the remaining budget for the current block row
                curr_per_block_row_budget = curr_per_block_row_budget - bits_used;
    
                % Append QP and encoded residues to the final output
                final_encodedResidues = [final_encodedResidues, QP, encodedResidues];
    
                % Store the quantized block for reconstruction
                quantizedResidues(row:block_row_end, col:block_col_end) = quantizedBlock;
    
                % Optionally, print the size of the encoded residues for debugging
                fprintf('Block at (%d, %d): QP=%d, Bits Used=%d, Remaining Budget=%d\n', ...
                        row, col, QP, bits_used, curr_per_block_row_budget);
            end
        end
        % Assign the final encoded residues to the output variable
        encodedResidues = final_encodedResidues;
    end

    % Prepend the frame type to the encoded data (as a header)
    encodedPredicitonModes = [frameTypeHeader, encodedPredicitonModes];
    encodedMotionVector = [frameTypeHeader, encodedMotionVector];
   
end


function reshaped_2d = reshape_3d_to_2d(matrix_3d)
    % Flatten the last two dimensions of the 3D matrix into a single dimension
    % Creating a 2D matrix where each row represents the flattened row blocks of the original 3D matrix
    [rows, cols, depth] = size(matrix_3d);
    reshaped_2d = reshape(matrix_3d, rows, cols * depth);
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


function decoded = rle_decode(encoded,datalength)
    % Decode data using a custom Run-Length Encoding (RLE) scheme
    % Return a 1D array containing the decoded values
    decoded = [];
    i = 1;
    
    while i <= length(encoded)
        if encoded(i) < 0
            % Negative value indicates a run of non-zero values
            non_zero_count = -encoded(i);
            i = i + 1;
            decoded = [decoded, encoded(i:i + non_zero_count - 1)];
            i = i + non_zero_count;
        elseif encoded(i) > 0
            % Positive value indicates a run of zeros
            zero_count = encoded(i);
            
            decoded = [decoded, zeros(1, zero_count)];
            i = i + 1;
        else
            % Zero value in encoded data indicates an actual zero in the original data
            decoded = [decoded, zeros(1, datalength-length(decoded))];
            i = length(encoded) + 1;
        end
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


function QP = findCorrectQP(per_block_row_budget, bitCountPerRow)
    % Iterate through the bitCountPerRow array
    for idx = 1:(length(bitCountPerRow))
        % Check if the budget falls within the current range
        if idx == length(bitCountPerRow)
            QP = idx - 1;
            return;
        end
        if per_block_row_budget >= bitCountPerRow(idx) && per_block_row_budget < bitCountPerRow(idx + 1)
            QP = idx - 1; % Return the corresponding QP value (adjust for 0-based indexing)
            return;
        end
    end
    % If no match is found, return -1 or an error message
    QP = -1; % Indicates that no valid QP was found
end

function Q = createQMatrix(blockSize, QP)
    if numel(blockSize) > 1
        rows = blockSize(1);
        cols = blockSize(2);
    else
        rows = blockSize;
        cols = blockSize;
    end
    
    Q = zeros(rows, cols);
    for x = 1:rows
        for y = 1:cols
            if (x + y <= rows)
                Q(x,y) = 2^QP;
            elseif (x + y == rows + 1)
                Q(x,y) = 2^(QP+1);
            else
                Q(x,y) = 2^(QP+2);
            end
        end
    end
end
