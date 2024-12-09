function data = exp_golomb_decode(encoded)
    % Decode a 1D array of binary values (0s and 1s) encoded using Exponential-Golomb coding
    data = [];  % Initialize an empty array to store decoded values
    idx = 1;
    
    while idx <= length(encoded)
        % Count leading zeros
        leading_zeros = 0;
        while idx <= length(encoded) && encoded(idx) == 0
            leading_zeros = leading_zeros + 1;
            idx = idx + 1;
        end
        
        % Read remainder of binary value
        if idx <= length(encoded) && encoded(idx) == 1
            idx = idx + 1;
            bin_code_length = leading_zeros;
            bin_code = [1, encoded(idx : idx + bin_code_length - 1)];
            idx = idx + bin_code_length;
            
            % Convert to decimal
            value = bin2dec(char(bin_code + '0'));
            
            % Map back to original value
            if mod(value, 2) == 0
                x = (value / 2);
            else
                x = -(value - 1) / 2;
            end
            
            data = [data, x];
        end
    end
end