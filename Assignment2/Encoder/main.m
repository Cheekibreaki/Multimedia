addpath('../Utils');  % For utils functions
addpath('../Outputs');  % For utils functions
addpath('../Decoder');  % For Decoder functions
% Parameters
filename = '../foreman_cif-1.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only_foreman.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y_foreman.yuv'; % File to store padded Y components
referenceFile = '../Outputs/referenceFrames.yuv';
decodedFile = '../Outputs/decoded_Y_foreman.yuv';
width = 352;                     % Frame width
height = 288;                    % Frame height
numFrames = 10;                 % Number of frames to process
searchRange = 4;                 % Search range r = 1,4, and 8
QP = 4;
j = 4;
VBSEnable = true;
FMEEnable = false;
FastME = false;

if(VBSEnable == true)
    j = j-1;
end

blockSize = 2^j;                   % Block size for motion estimation
dct_blockSize = 2^j;
I_Period = 8; 
nRefFrames = 4;                 % Can take value from 1 to 4

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

% % encoder
encoder(referenceFile, paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable, FMEEnable,FastME );
[total_byte,bytes_list] = decoder(decodedFile);
% %decoder
% compareYUVFrames(referenceFile, outputFile, decodedFile, width, height, numFrames);
% 
% 
% 
% calculatePSNR(decodedFile, paddedOutputFile, width, height, numFrames)



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

