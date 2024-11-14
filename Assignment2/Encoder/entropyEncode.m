function [encodedMotionVector,encodedPredicitonModes,encodedResidues] = entropyEncode(frame_type, motionVector3d, predicitonModes2d, residues2d,vbs_matrix)
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
                    % [motion_block,previous_motion_vector_block] = diffEncoding_block(motion_block,'mv',previous_motion_vector_block);
                    
                    if all(vbs_block(:) == 0)  % If this 2x2 block in vbs_matrix is all zeros
                        % Fetch the entire 2x2x3 block from motionVector3d
                        motionVector2d = reshape_3d_to_2d(motion_block);
                        motionVector1d = zigzag(motionVector2d);
                        motionVector1d = [0, motionVector1d];  % Prefix with 0 to indicate all zeros
                        resultmotionVector1d = [resultmotionVector1d, motionVector1d];  % Append to the result
                    elseif all(vbs_block(:) == 1)  % If this 2x2 block in vbs_matrix is all ones
                        % Fetch the 1x1x3 top-left element of the 2x2x3 motionVector3d block
                        top_left_block = motion_block(1, 1, :);
                        motionVector2d = reshape_3d_to_2d(top_left_block);
                        motionVector1d = zigzag(motionVector2d);
                        motionVector1d = [1, motionVector1d];  % Prefix with 1 to indicate all ones
                        resultmotionVector1d = [resultmotionVector1d, motionVector1d];  % Append to the result
                    else
                        error('ErrorID:1', 'a block can not contain both 1 and 0!');
                        %This shouldn't happen! something goes wrong
                    end
                end
            end
    
            % Encode the result using Exp-Golomb encoding
            encodedMotionVector = exp_golomb_encode(resultmotionVector1d);
            a=1
        else
            % For block encoding if vbs_matrix does not exist
            motionVector2d = reshape_3d_to_2d(motionVector3d);
            motionVector1d = zigzag(motionVector2d);
            encodedMotionVector = exp_golomb_encode(motionVector1d);
        end

    elseif frame_type == 1
        predicitonModes1d = zigzag(predicitonModes2d);
        encodedPredicitonModes = exp_golomb_encode(predicitonModes1d); 

    end

        residues1d = zigzag(residues2d);
        residuesRLE = rle_encode(residues1d); 
        encodedResidues = exp_golomb_encode(residuesRLE);

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

% function encoded = exp_golomb_encode(data)
%     % Encode data using Exponential-Golomb coding based on specific rules
%     % and return a 1D cell array of binary strings
%     encoded = cell(1, length(data));  % Initialize a cell array to store binary strings
%     for i = 1:length(data)
%         x = data(i);
%         if x > 0
%             value = 2 * x;  % Positive x to odd integer (2x-1)
%         else
%             value = (-2 * x)+1;     % Non-positive x to even integer (-2x)
%         end
%         bin_code = dec2bin(value);  % Convert to binary
%         leading_zeros = floor(log2(value)) + 1;  % Compute the number of leading zeros
%         code = [repmat('0', 1, leading_zeros - 1), '1', bin_code(2:end)];  % Assemble the code
%         encoded{i} = code;  % Store binary code in cell array
%     end
% end


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


% testdata = [-31, 9, -4, 8, 1, -3, 4, 4, 2, 4, 0, 4, 0, 0, -4, 0, 0, 1, 0, 0]
% encoded = rle_encode(testdata)
% decoded = rle_decode(encoded,20)