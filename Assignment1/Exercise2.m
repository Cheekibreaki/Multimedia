function processYUV(filename, width, height, numFrames)
    % Parameters:
    % filename - path to the .yuv file
    % width - width of the video
    % height - height of the video
    % numFrames - number of frames to process

    % Open the YUV file
    fid = fopen(filename, 'r');

    for frameIdx = 1:numFrames
        % Read Y plane
        Y = fread(fid, [width, height], 'uint8')';

        % Dump Y-only frame to file
        yFilename = sprintf('Y_only_frame_%d.yuv', frameIdx);
        fidY = fopen(yFilename, 'w');
        fwrite(fidY, Y, 'uint8');
        fclose(fidY);
        
        % Process with block sizes 2x2, 8x8, 64x64
        blockSizes = [8];
        for blockSize = blockSizes
            processYBlock(Y, blockSize, width, height, frameIdx);
        end
    end

    % Close the original file
    fclose(fid);
end

function processYBlock(Y, blockSize, width, height, frameIdx)
    % Parameters:
    % Y - Y component frame
    % blockSize - size of the block (2, 8, 64)
    % width - width of the frame
    % height - height of the frame
    % frameIdx - current frame number

    % Pad the frame if necessary
    paddedY = padFrame(Y, blockSize, width, height);

    % Block processing
    [paddedHeight, paddedWidth] = size(paddedY);
    newY = paddedY;

    for row = 1:blockSize:paddedHeight
        for col = 1:blockSize:paddedWidth
            % Extract the block
            block = paddedY(row:min(row+blockSize-1, paddedHeight), col:min(col+blockSize-1, paddedWidth));

            % Calculate the average value of the block
            avgValue = round(mean(block(:)));

            % Replace the block with the average value
            newY(row:min(row+blockSize-1, paddedHeight), col:min(col+blockSize-1, paddedWidth)) = avgValue;
        end
    end

    % Save the Y-only block-averaged frame to file
    yBlockFilename = sprintf('Y_only_block_averaged_%dx%d_frame_%d.yuv', blockSize, blockSize, frameIdx);
    fidBlock = fopen(yBlockFilename, 'w');
    fwrite(fidBlock, newY', 'uint8');  % Transpose for correct orientation
    fclose(fidBlock);

    % Display subjective comparison (original vs. block averaged)
    figure;
    subplot(1, 2, 1);
    imshow(uint8(Y));
    title(['Original Y Frame ', num2str(frameIdx)]);
    
    subplot(1, 2, 2);
    imshow(uint8(newY));
    title(['Block Averaged Y Frame ', num2str(frameIdx), ' (', num2str(blockSize), 'x', num2str(blockSize), ')']);
    
    % Frame differencing and comparison
    difference = abs(double(paddedY) - double(newY));
    differenceMagnified = difference * 10;  % Magnify the deltas
    figure;
    imshow(uint8(differenceMagnified));
    title(['Difference Magnified for Frame ', num2str(frameIdx), ' (', num2str(blockSize), 'x', num2str(blockSize), ')']);
    
    % Compute PSNR and SSIM
    psnrValue = psnr(uint8(newY), uint8(paddedY));
    ssimValue = ssim(uint8(newY), uint8(paddedY));
    
    fprintf('PSNR (Block size %dx%d): %.2f dB\n', blockSize, blockSize, psnrValue);
    fprintf('SSIM (Block size %dx%d): %.4f\n', blockSize, blockSize, ssimValue);
end

function paddedY = padFrame(Y, blockSize, width, height)
    % Pad the frame with gray (128) if necessary
    padRight = mod(width, blockSize);
    padBottom = mod(height, blockSize);

    if padRight ~= 0
        padRight = blockSize - padRight;
    end
    if padBottom ~= 0
        padBottom = blockSize - padBottom;
    end

    paddedY = padarray(Y, [padBottom, padRight], 128, 'post');
end

function plotResults(psnrValues, ssimValues, blockSizes)
    % Plot PSNR and SSIM vs block sizes
    figure;
    subplot(1, 2, 1);
    plot(blockSizes, psnrValues, '-o');
    xlabel('Block Size');
    ylabel('Average PSNR (dB)');
    title('PSNR vs Block Size');
    
    subplot(1, 2, 2);
    plot(blockSizes, ssimValues, '-o');
    xlabel('Block Size');
    ylabel('Average SSIM');
    title('SSIM vs Block Size');
end

processYUV('foreman_cif-1.yuv',352,288,2)