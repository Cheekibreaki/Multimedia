function compareYUVFrames(referenceFile, originalFile, decodedFile, width, height, numFrames)
    % compareYUVFrames: This function compares frames from the original
    % YUV file and the decoded YUV file by displaying them side-by-side.
    %
    % Parameters:
    %   originalFile - Path to the original YUV file
    %   decodedFile  - Path to the decoded YUV file
    %   width        - Width of the frames
    %   height       - Height of the frames
    %   numFrames    - Number of frames to compare

    % Open the original and decoded files for reading
    fidReference = fopen(referenceFile, 'r');
    fidDecoded = fopen(decodedFile, 'r');
    
    fread(fidReference, [width, height], 'uint8')';

    % Iterate through each frame to compare
    for frameIdx = 1:numFrames
        % Read the Y component of the original frame
        
        refFrame = fread(fidReference, [width, height], 'uint8')';
        
        % Read the Y component of the decoded frame
        decodedFrame = fread(fidDecoded, [width, height], 'uint8')';
        
        % Calculate the absolute difference between the frames
        differenceFrame = abs(double(refFrame) - double(decodedFrame));

        % Display the frames for comparison
        figure;
        
        % Display the original frame
        subplot(1, 3, 1);
        imshow(refFrame, []);
        title(sprintf('Reference Frame %d', frameIdx));
        
        % Display the decoded frame
        subplot(1, 3, 2);
        imshow(decodedFrame, []);
        title(sprintf('Decoded Frame %d', frameIdx));
        
        % Display the difference frame
        subplot(1, 3, 3);
        imshow(differenceFrame, []);
        title(sprintf('Difference Frame %d', frameIdx));
        
        % Pause to allow viewing before moving to the next frame
        pause(1);  % Adjust or remove pause as needed
    end
    
    % Close the files
    fclose(fidReference);
    fclose(fidDecoded);
end