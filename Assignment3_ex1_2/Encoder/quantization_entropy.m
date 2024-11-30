function [quantizedResiduals] = quantization_entropy(residuals, dct_blockSize, width, height, baseQP,RCflag,per_block_row_budget, bitCountPerRow, vbs_matrix)
    quantizedResiduals = zeros(size(residuals));
    if exist('vbs_matrix', 'var') && ~isempty(vbs_matrix)
    [vbsRows, vbsCols] = size(vbs_matrix);
    


    %some thing is not right with baseqp adjustment









    row_bits_used = per_block_row_budget;
    % Iterate over the vbs_matrix in steps of 2 to process 2x2 blocks
    for blockY = 1:2:vbsRows
        if RCflag
            next_row_budget = per_block_row_budget + (per_block_row_budget - row_bits_used);
            baseQP = findCorrectQP(next_row_budget,bitCountPerRow);
            row_bits_used = 0;
        end
        for blockX = 1:2:vbsCols
            % Extract the 2x2 block from vbs_matrix
            vbs_block = vbs_matrix(blockY:blockY+1, blockX:blockX+1);

            % Determine the starting position in the residuals
            rowOffset = (blockY - 1) * dct_blockSize + 1;
            colOffset = (blockX - 1) * dct_blockSize + 1;

            % Calculate the actual block size, accounting for image boundaries
            actualBlockHeight = min(2 * dct_blockSize, height - rowOffset + 1);
            actualBlockWidth = min(2 * dct_blockSize, width - colOffset + 1);

            % Check if the vbs_block is all ones or contains zeros
            if all(vbs_block(:) == 0)
                % Process as a large block
                % Extract the residual block
                block = residuals(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1);

                % Apply DCT to the block
                dctBlock = dct2(double(block));

                % Create the quantization matrix
                Q = createQMatrix(size(dctBlock), baseQP);

                % Quantize the DCT coefficients
                quantizedBlock = round(dctBlock ./ Q);

                % Store the quantized block
                quantizedResiduals(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = quantizedBlock;
            else
                % Process each sub-block separately
                for subBlockY = 0:1
                    for subBlockX = 0:1
                        % Calculate sub-block offsets
                        subRowOffset = rowOffset + subBlockY * dct_blockSize;
                        subColOffset = colOffset + subBlockX * dct_blockSize;

                        % Check if within image boundaries
                        if subRowOffset > height || subColOffset > width
                            continue;
                        end

                        % Calculate actual sub-block size
                        actualSubBlockHeight = min(dct_blockSize, height - subRowOffset + 1);
                        actualSubBlockWidth = min(dct_blockSize, width - subColOffset + 1);

                        % Extract the residual sub-block
                        block = residuals(subRowOffset:subRowOffset+actualSubBlockHeight-1, subColOffset:subColOffset+actualSubBlockWidth-1);

                        % Apply DCT to the sub-block
                        dctBlock = dct2(double(block));

                        % Create the quantization matrix
                        Q = createQMatrix(size(dctBlock), baseQP -1);

                        % Quantize the DCT coefficients
                        quantizedBlock = round(dctBlock ./ Q);

                        % Store the quantized sub-block
                        quantizedResiduals(subRowOffset:subRowOffset+actualSubBlockHeight-1, subColOffset:subColOffset+actualSubBlockWidth-1) = quantizedBlock;

                    end
                end
            end
            if RCflag
                 % Flatten the quantized block using zigzag scan
                quantizedBlock1d = zigzag(quantizedBlock);
                residuesRLE = rle_encode(quantizedBlock1d); 
                encodedResidues = exp_golomb_encode(residuesRLE);
                encodedResidues_length = length(encodedResidues);
                % Calculate bits used by this block
                row_bits_used = row_bits_used + encodedResidues_length;
            end
        end
    end
    else
        quantizedResiduals = zeros(size(residuals));
        for row = 1:dct_blockSize:height
            for col = 1:dct_blockSize:width
                block = residuals(row:row+dct_blockSize-1, col:col+dct_blockSize-1);
                dctBlock = dct2(double(block));

                Q = createQMatrix(size(block), baseQP);
                quantizedBlock = round(dctBlock ./ Q);
                quantizedResiduals(row:row+dct_blockSize-1, col:col+dct_blockSize-1) = quantizedBlock;
            end
            if RCflag
                 % Flatten the quantized block using zigzag scan
                quantizedBlock1d = zigzag(quantizedBlock);
                residuesRLE = rle_encode(quantizedBlock1d); 
                encodedResidues = exp_golomb_encode(residuesRLE);
    
                % Calculate bits used by this block
                row_bits_used = bits_used + length(encodedResidues);
            end
        end
    end
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
    % Convert variables to int32 if necessary
    per_block_row_budget = int32(per_block_row_budget);
    bitCountPerRow = int32(bitCountPerRow);
    
    % Iterate through the bitCountPerRow array up to the second last element
    for idx = 1:(length(bitCountPerRow) - 1)
        val = bitCountPerRow(idx + 1);
        next_val = bitCountPerRow(idx);
        % Check if the budget falls within the current range
        if (per_block_row_budget > val && per_block_row_budget <= next_val)
            QP = idx - 1; % Adjust for 0-based indexing
            return;
        end
    end
    % Check if per_block_row_budget is greater than or equal to the last element
    if per_block_row_budget >= bitCountPerRow(end)
        QP = length(bitCountPerRow) - 1; % Adjust for 0-based indexing
        return;
    end
    % If no match is found, return -1
    if(per_block_row_budget > bitCountPerRow(1))
        QP = 0;
    else
        QP = 9;
    end
end


