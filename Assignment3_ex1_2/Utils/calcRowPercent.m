function total_per_row_perc = calculateRowPercent(total_per_row_bits_used, total_bits_used)
    % Function to calculate the percentage of bits used per row
    %
    % Inputs:
    %   - total_per_row_bits_used: Array of bits used per row
    %   - total_bits_used: Total number of bits used
    %
    % Outputs:
    %   - total_per_row_perc: Array of percentages for each row
    
    % Validate inputs
    if total_bits_used == 0
        total_per_row_perc = [];
        return;
    end

    % Calculate percentages (vectorized for performance)
    total_per_row_perc = total_per_row_bits_used / total_bits_used;
end