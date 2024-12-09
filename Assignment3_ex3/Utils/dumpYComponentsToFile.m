function dumpYComponentsToFile(filename, width, height, numFrames, outputFile)
    % Parameters:
    % filename - path to the .yuv file
    % width - width of the video
    % height - height of the video
    % numFrames - number of frames to process
    % outputFile - path to the file where Y components will be dumped
    
    % Open the YUV file
    fid = fopen(filename, 'r');
    uvSize = width * height /4;
    
    % Open the output file for writing the Y components
    fidOut = fopen(outputFile, 'w');
    
    for frameIdx = 1:numFrames
        % Read Y plane from the YUV file
        Y = fread(fid, [width, height], 'uint8');
        
        % Dump the Y-component of the frame into the output file
        fwrite(fidOut, Y, 'uint8');

        %skip the U and V components
        fseek(fid, 2 * uvSize,'cof');
    end
    
    % Close both the input YUV file and the output file
    fclose(fid);
    fclose(fidOut);
    
    fprintf('Dumped Y components of %d frames into file: %s\n', numFrames, outputFile);
end