addpath('../Utils');  % For utils functions
addpath('../Outputs');  % For utils functions
addpath('../Decoder');  % For Decoder functions
% Parameters
filename = '../CIF.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y.yuv'; % File to store padded Y components
referenceFile = '../Outputs/referenceFrames.yuv';
decodedFile = '../Outputs/decoded_Y.yuv';
predictedFile = '../Outputs/predictedFrames.yuv';
width = 352;                     % Frame width
height = 288;                    % Frame height
numFrames = 10;                 % Number of frames to process
searchRange = 4;                 % Search range r = 1,4, and 8
QP = 4;
j = 4;
VBSEnable = true;
FMEEnable = true;
FastME = true;
mode =2;

if(VBSEnable == true)
    j = j-1;
end

blockSize = 2^j;                   % Block size for motion estimation
dct_blockSize = 2^j;
I_Period = 8; 
nRefFrames = 1;                 % Can take value from 1 to 4

function lambda = get_lambda_for_qp(QP)
   if QP == 1
            lambda = 0.02;
        elseif QP == 2
            lambda = 0.1667;
        elseif QP == 4
            lambda = 0.1667;
        elseif QP == 7
            lambda = 0.222;
        elseif QP == 10
            lambda = 0.5;
        else
            % Handle other cases if needed
            lambda = NaN;  % Example: set to NaN if QP is not 1, 4, 7, or 10
   end
end
lambda = get_lambda_for_qp(QP);
% Pre-process

QPs = [1,2,4,7,10];


dumpYComponentsToFile(filename, width, height, numFrames, outputFile);

[paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

coreCount = 2;
% Start parallel pool
parpool(coreCount);
    
%  encoder
if mode == 1
    encoder_mode1(referenceFile, paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable, FMEEnable,FastME,mode);
elseif mode == 3
    encoder_mode3(referenceFile, paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable, FMEEnable,FastME,mode);
else
    encoder(referenceFile, paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable, FMEEnable,FastME,mode);
end
[total_byte,bytes_list] = decoder(decodedFile,mode);
% decoder
 compareYUVFrames(predictedFile, outputFile, decodedFile, width, height, numFrames);
% 
% 
% 
% calculatePSNR(decodedFile, paddedOutputFile, width, height, numFrames)

% Delete the parallel pool
    delete(gcp('nocreate'));

%All the graph for report:
%generate_rd_analysis();

% 
% 
% analyze_vbs_statistics(filename, outputFile, paddedOutputFile, referenceFile, decodedFile, ...
%                        width, height, numFrames, blockSize, searchRange, ...
%                        I_Period, lambda, VBSEnable, FMEEnable, FastME, ...
%                        nRefFrames, QPs);
% 
% lambdatest(filename, outputFile, paddedOutputFile, referenceFile, decodedFile, ...
%                        width, height, numFrames, blockSize, searchRange, ...
%                        I_Period, lambda, VBSEnable, FMEEnable, FastME, ...
%                        nRefFrames, QPs);

%analyze_reference_frames();

%analyze_lambda_rd_plot();
