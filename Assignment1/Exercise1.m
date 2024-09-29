function upscaleYUV420to444(filename, width, height, numFrames, isNoise)
    % Parameters:
    % filename - path to the .yuv file
    % width - width of the video
    % height - height of the video
    % numFrames - number of frames to process
    % isNoise - specifies which component to add noise to: 'Y', 'U', 'V', 'R', 'G', 'B'

    % Calculate sizes for Y, U, and V planes
    ySize = width * height;
    uvSize = ySize / 4;  % U and V planes are quarter the size of Y in 4:2:0 format

    % Open the YUV file
    fid = fopen(filename, 'r');

    for frameIdx = 1:numFrames
        % Read Y plane
        Y = fread(fid, [width, height], 'uint8');

        % Read U and V planes (quarter resolution)
        U = fread(fid, [width/2, height/2], 'uint8');
        V = fread(fid, [width/2, height/2], 'uint8');

        % Upscale U and V planes to full resolution (bilinear interpolation)
        U_full = imresize(U, [width, height], 'bilinear');
        V_full = imresize(V, [width, height], 'bilinear');

        % Add random noise to Y, U, or V planes if specified
        if strcmp(isNoise, 'Y')
            Y = addRandomNoise(Y);
        elseif strcmp(isNoise, 'U')
            U_full = addRandomNoise(U_full);
        elseif strcmp(isNoise, 'V')
            V_full = addRandomNoise(V_full);
        end

        % Combine Y, U_full, and V_full into a YUV 4:4:4 frame
        YUV_444 = cat(3, Y, U_full, V_full);

        % Convert to RGB (use the modified yuv2rgb function)
        RGB = yuv2rgb(YUV_444, isNoise);

        % Display the frame
        imshow(RGB);
        title(['Frame ', num2str(frameIdx)]);
        pause(0.1);  % Pause to view each frame
    end

    % Close the file
    fclose(fid);
end

function RGB = yuv2rgb(YUV, isNoise)
    % Convert YUV to RGB using standard formula
    % YUV is assumed to be in the range 0-255
    % isNoise - specifies which component to add noise to: 'R', 'G', 'B'

    Y = double(YUV(:,:,1))';
    U = double(YUV(:,:,2))';
    V = double(YUV(:,:,3))';

    % Standard conversion matrix from YUV to RGB
    R = 1.164*(Y-16) + 1.596*(V-128);
    G = 1.164*(Y-16) - 0.392*(U-128) - 0.813*(V-128);
    B = 1.164*(Y-16) + 2.017*(U-128);

    % Clip values to [0, 255]
    R = max(min(R, 255), 0);
    G = max(min(G, 255), 0);
    B = max(min(B, 255), 0);

    % Add random noise to R, G, or B components if specified
    if strcmp(isNoise, 'R')
        R = addRandomNoise(R);
    elseif strcmp(isNoise, 'G')
        G = addRandomNoise(G);
    elseif strcmp(isNoise, 'B')
        B = addRandomNoise(B);
    end

    % Combine into RGB image
    RGB = cat(3, uint8(R), uint8(G), uint8(B));
end

function componentWithNoise = addRandomNoise(component)
    % Add random noise to the component
    noiseLevel = 10;  % You can adjust the noise level as needed
    noise = noiseLevel * randn(size(component));  % Gaussian noise
    componentWithNoise = double(component) + noise;
    
    % Clip values to [0, 255] to keep within valid range
    componentWithNoise = max(min(componentWithNoise, 255), 0);
end


upscaleYUV420to444('foreman_cif-1.yuv',352,288,100,'Y')

upscaleYUV420to444('foreman_cif-1.yuv',352,288,100,'U')
upscaleYUV420to444('foreman_cif-1.yuv',352,288,100,'V')
upscaleYUV420to444('foreman_cif-1.yuv',352,288,100,'R')
upscaleYUV420to444('foreman_cif-1.yuv',352,288,100,'G')
upscaleYUV420to444('foreman_cif-1.yuv',352,288,100,'B')

