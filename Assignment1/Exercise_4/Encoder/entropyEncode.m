function [encodedMotionVector,encodedPredicitonModes,encodedResidues] = entropyEncode(frame_type, motionVector3d, predicitonModes2d, residues2d)
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

        motionVector2d = reshape_3d_to_2d(motionVector3d);
        motionVector1d = zigzag(motionVector2d);
        motionVectorRLE = rle_encode(motionVector1d);
        encodedMotionVector = exp_golomb_encode(motionVectorRLE);

    elseif frame_type == 1
        predicitonModes1d = zigzag(predicitonModes2d);
        predicitonModesRLE = rle_encode(predicitonModes1d);
        encodedPredicitonModes = exp_golomb_encode(predicitonModesRLE); 

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