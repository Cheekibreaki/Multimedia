function total_per_row_qp = findQP(total_per_row_perc, baseQP)
    % findQP - Determines QP values for the second pass based on the first pass statistics.
    %
    % Inputs:
    %   total_per_row_perc: A vector of normalized bit percentages per row from the first pass.
    %   baseQP: The base QP used in the first pass.
    %
    % Output:
    %   total_per_row_qp: A vector of integer QP values for each row during the second pass.

    % Compute the baseline
    baseline = 1 / length(total_per_row_perc);

    % Initialize output array
    total_per_row_qp = [];

    % Define a scaling factor
    scaling_factor = 5.0;

    for row_perc = total_per_row_perc
        if row_perc >= baseline
            % Increase QP for rows that exceeded the baseline bit usage
            newQP = baseQP + (row_perc - baseline) * scaling_factor * baseQP;
        else
            % Decrease QP for rows that used fewer bits than baseline
            newQP = baseQP - (baseline - row_perc) * scaling_factor * baseQP;
            newQP = max(newQP, 0); % Ensure QP does not go below zero
        end

        % Round to nearest integer
        newQP = round(newQP);

        % Append to output array
        total_per_row_qp = [total_per_row_qp; newQP];
    end
end
