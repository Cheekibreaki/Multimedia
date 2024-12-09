function interpolatedFrame = interpolateFrame(referenceFrame)
    % creating a 2X resolution frame using bilinear interpolation.
    %
    % Parameters:
    %   referenceFrame - The original reference frame (grayscale image)
    %
    % Returns:
    %   interpolatedFrame - The frame with interpolated values (2X resolution)

    [height, width] = size(referenceFrame);

    % Initialize the interpolated frame with 2X the resolution
    interpolatedFrame = zeros(2 * height -1 , 2 * width -1);

    % Fill in the actual pixels at even indices
    interpolatedFrame(1:2:end, 1:2:end) = referenceFrame;

    % Interpolate horizontally
    for y = 1:2:2 * height -1
        for x = 2:2:2 * width - 1
            interpolatedFrame(y, x) = round((interpolatedFrame(y, x - 1) + interpolatedFrame(y, x + 1)) / 2);
        end
    end

    % Interpolate vertically
    for y = 2:2:2 * height - 1
        for x = 1:2 * width - 1
            interpolatedFrame(y, x) = round((interpolatedFrame(y - 1, x) + interpolatedFrame(y + 1, x)) / 2);
        end
    end

    % Interpolate diagonally
    for y = 2:2:2 * height - 1
        for x = 2:2:2 * width - 1
            interpolatedFrame(y, x) = round((interpolatedFrame(y - 1, x - 1) + interpolatedFrame(y - 1, x + 1) + ...
                                       interpolatedFrame(y + 1, x - 1) + interpolatedFrame(y + 1, x + 1)) / 4);
        end
    end
end
