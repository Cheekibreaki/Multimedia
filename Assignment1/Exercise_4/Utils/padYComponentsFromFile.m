function [paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile)
    % Parameters:
    % outputFile - path to the Y-only file that contains all Y components
    % numFrames - number of Y-only frames to process
    % width - width of each Y frame
    % height - height of each Y frame
    % blockSizes - array of block sizes (e.g., [2, 8, 64])
    
    % Open the Y-only file
    fidY = fopen(outputFile, 'r');

    % Open a file to save the padded Y frames
    fidPadded = fopen(paddedOutputFile, 'w');
    
    for frameIdx = 1:numFrames
        % Read the Y component of the current frame
        Y = fread(fidY, [width, height], 'uint8')';
       
        [paddedY,paddedWidth,paddedHeight] = padFrame(Y, width, height, blockSize);
       
        fwrite(fidPadded, paddedY', 'uint8');                
        
    end
    
    % Close the files after processing all frames
    fclose(fidY);
    fclose(fidPadded);
    fprintf('Padded Y components of %d frames into file: %s\n', numFrames, paddedOutputFile);
end



function [paddedY,paddedWidth,paddedHeight] = padFrame(Y, width, height, blockSize)
    % Pad the frame with gray (128) if necessary to make it divisible by blockSize
    padRight = mod(width, blockSize);
    padBottom = mod(height, blockSize);

    if padRight ~= 0
        padRight = blockSize - padRight;
    end
    if padBottom ~= 0
        padBottom = blockSize - padBottom;
    end
    
    paddedWidth = width + padRight;
    paddedHeight = height + padBottom;

    % Pad the frame on the right and bottom with 128 (gray)
    paddedY = padarray(Y, [padRight, padBottom], 128, 'post');

end

