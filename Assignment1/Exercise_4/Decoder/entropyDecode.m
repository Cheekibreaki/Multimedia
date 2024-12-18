function [decodedMotionVector3d,decodedPredicitonModes2d,decodedResidues2d] = entropyDecode(frame_type, encodedMotionVector1d, encodedPredicitonModes1d, encodedResidues1d,mvwidth, mvheight,  predwidth, predheight,  reswidth, resheight)
    
    %switch the order of width & height
    temp = mvwidth;
    mvwidth = mvheight;
    mvheight = temp;

    temp = predwidth;
    predwidth = predheight;
    predheight = temp;
    
    temp = reswidth;
    reswidth = resheight;
    resheight = temp;

    % Input:
    % frame_type: 1 for I-frame, 0 for P-frame 
    % pred_diff: array containing differential prediction information (modes for intra or motion vectors for inter)
    % dct_coeffs: block of quantized DCT coefficients
    % block_size: size of the block (e.g., 4 for 4x4 blocks)

    decodedMotionVector3d = [];
    decodedPredicitonModes2d = [];
    decodedResidues2d = [];
    if frame_type == 0  % Assuming data handling for P-frames
        
        motionVectorRevEGC = exp_golomb_decode(encodedMotionVector1d);
        motionVector2d = invzigzag(motionVectorRevEGC, mvheight, length(motionVectorRevEGC)/mvheight);
        decodedMotionVector3d = reshape_2d_to_3d(motionVector2d, mvwidth, mvheight, 2);

    elseif frame_type == 1
    
        predictedModeRevEGC = exp_golomb_decode(encodedPredicitonModes1d); 
        decodedPredicitonModes2d = invzigzag(predictedModeRevEGC, predwidth, predheight);
        
    end
        
        residuesRevEGC = exp_golomb_decode(encodedResidues1d);
        decodedResidues1d = rle_decode(residuesRevEGC, reswidth * resheight); 
        decodedResidues2d = invzigzag(decodedResidues1d, reswidth, resheight);
    
end


function matrix_3d = reshape_2d_to_3d(reshaped_2d, rows, cols, depth)
    % Reshape the 2D matrix back to the original 3D matrix
    % Each row block in the original 3D matrix is represented in flattened form in reshaped_2d
    matrix_3d = reshape(reshaped_2d, cols, rows, depth);
end

function matrix = invzigzag(zigzag_order, rows, cols)
    % Inverse Zigzag function to reorder a 1D array back into a 2D matrix
    matrix = zeros(rows, cols);
    index = 1;
    for s = 1:(rows + cols - 1)
        % Even diagonals: top-right to bottom-left
        for i = max(1, s - cols + 1):min(rows, s)
            j = s - i + 1;
            matrix(i, j) = zigzag_order(index);
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

function data = exp_golomb_decode(encoded)
    % Decode a 1D array of binary values (0s and 1s) encoded using Exponential-Golomb coding
    data = [];  % Initialize an empty array to store the decoded values
    idx = 1;
    
    while idx <= length(encoded)
        % Count the number of leading zeros
        leading_zeros = 0;
        while idx <= length(encoded) && encoded(idx) == 0
            leading_zeros = leading_zeros + 1;
            idx = idx + 1;
        end
        
        % Read the remainder of the binary value
        if idx <= length(encoded) && encoded(idx) == 1
            idx = idx + 1;
            bin_code_length = leading_zeros;
            bin_code = [1, encoded(idx : idx + bin_code_length - 1)];
            idx = idx + bin_code_length;
            
            % Convert binary code to decimal value
            value = bin2dec(char(bin_code + '0'));
            
            % Map value back to original data
            if mod(value, 2) == 0
                x = (value / 2);
            else
                x = - (value - 1) / 2;
            end
            
            data = [data, x];  % Append the decoded value to the output array
        end
    end
end


function decoded = rle_decode(encoded,datalength)
    % Decode data using a custom Run-Length Encoding (RLE) scheme
    % Return a 1D array containing the decoded values of shape 1x(blockSize * blockSize)
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