function avgPsnr = calculatePSNR(decodedFile, originalFile, width, height, numFrames)
    % Parameters:
    %   originalFile - Path to the original Y only file
    %   decodedFile  - Path to the decoded Y only file
    %   width        - Width of the frames
    %   height       - Height of the frames
    %   numFrames    - Number of frames to compare

    % Open the original and decoded files for reading
    fidOriginal = fopen(originalFile, 'r');
    fidDecoded = fopen(decodedFile, 'r');
    
    totalPsnr = 0;
    
    % Iterate through each frame to compare
    for frameIdx = 1:numFrames
        % Read frames
        originalFrame = fread(fidOriginal, [width, height], 'uint8');
        decodedFrame = fread(fidDecoded, [width, height], 'uint8');
        

        % Calculate MSE
        mse = mean((double(originalFrame(:)) - double(decodedFrame(:))).^2);

        % Calculate PSNR
        if mse > 0
            psnrValue = 10 * log10(255^2 / mse);
        else
            psnrValue = Inf; % Perfect match
        end
        % psnrValue = psnr(originalFrame,decodedFrame)
        
        % Accumulate PSNR
        totalPsnr = totalPsnr + psnrValue;
    end
    
    % Average PSNR over all frames
    avgPsnr = totalPsnr / numFrames;
    %fprintf('Average PSNR: %.2f dB\n', avgPsnr);
    
    % Close the files
    fclose(fidOriginal);
    fclose(fidDecoded);
end