function lambdatest(filename, outputFile, paddedOutputFile, referenceFile, decodedFile, ...
                       width, height, numFrames, blockSize, searchRange, ...
                       I_Period, lambda, VBSEnable, FMEEnable, FastME, ...
                       nRefFrames, QPs)
    dct_blockSize = blockSize
    meanPSNRValues = [];
    for QP = QPs
        VBSEnable = true;
        
        % https://ieeexplore.ieee.org/document/1626308
        if QP == 1
            Lambdas = linspace(0, 1, 7);
        elseif QP == 2
            Lambdas = linspace(0, 1, 7);
        elseif QP == 4
            Lambdas = linspace(0, 1, 7);
        elseif QP == 7
            Lambdas = linspace(0, 1, 10);
        elseif QP == 10
            Lambdas = linspace(0, 1, 7);
        end
        

       % λMODE=0.85×2(QP−12)/3
        i = 1

        for Lambda = Lambdas


            dumpYComponentsToFile(filename, width, height, numFrames, outputFile);

            [paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);
            
            % encoder
            encoder(referenceFile, paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,Lambda,VBSEnable, FMEEnable,FastME );
            %decoder
            [total_byte,bytes_list] = decoder(decodedFile);
            
            
            psnrValues(i) = calculatePSNR(decodedFile, paddedOutputFile, width, height, numFrames);
            total_bits(i) = total_byte * 8;
           
           
            % disp(strcat('Lambda=', string(Lambda), ', PSNR=', string(meanPsnrValue), ', total Bits=', string(totalBitsFrames), ', PSNR/totalBits=', string(meanPsnrValue/meanTotal_bits)));
            i = i + 1
        

        
        
        
        
        end

         plotAgainstFrame(total_bits, psnrValues, 'totalBits', 'PSNR', strcat('R-D plot when QP=',string(QP), ' varying Lambda'),Lambdas);
    end
end




function plotAgainstFrame(x_frame, yVals, xaxisLabel, yaxisLabel, titleStr, lambdas)
    figure;

    % Plot the data
    for i = 1:size(yVals, 1)
        plot(x_frame, yVals(i, :), '-o', 'DisplayName', strcat('QP-', num2str(i)));
        hold on;
    end

    xlabel(xaxisLabel);
    ylabel(yaxisLabel);
    grid on;
    hold off;
    title(titleStr);

    % Add a legend
    legend('show');

    % Customize data tips
    dcm = datacursormode(gcf);
    set(dcm, 'UpdateFcn', @(obj, event) dataTipCallback(obj, event, x_frame, yVals, lambdas));
end

function txt = dataTipCallback(~, event, x_frame, yVals, lambdas)
    % Get the data cursor position
    pos = get(event, 'Position');
    x = pos(1);
    y = pos(2);

    % Find the closest point in x_frame
    [~, idx] = min(abs(x_frame - x));

    % Retrieve the corresponding lambda value
    lambdaValue = lambdas(idx);

    % Construct the text to display
    txt = {['X: ', num2str(x)], ...
           ['Y: ', num2str(y)], ...
           ['Lambda: ', num2str(lambdaValue)]};
end