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
numFrames = 3;                 % Number of frames to process
searchRange = 8;                 % Search range r = 1,4, and 8
QP = 3;
j = 4;
VBSEnable = true;
FMEEnable = false;
FastME = false;

if(VBSEnable == true)
    j = j-1;
end

blockSize = 2^j;                   % Block size for motion estimation
dct_blockSize = 2^j;
I_Period = 5; 
nRefFrames = 4;                 % Can take value from 1 to 4

lambda = 0.65;
% Pre-process

QPs = [2,4,7];


% dumpYComponentsToFile(filename, width, height, numFrames, outputFile);
% 
% [paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);
% 
% % encoder
% encoder(referenceFile, paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable, FMEEnable,FastME );
% [total_byte,bytes_list] = decoder(decodedFile);
% %decoder
% compareYUVFrames(referenceFile, outputFile, decodedFile, width, height, numFrames);
% 
% 
% 
% calculatePSNR(decodedFile, paddedOutputFile, width, height, numFrames)



%All the graph for report:
%generate_rd_analysis();

% 
% analyze_vbs_statistics(filename, outputFile, paddedOutputFile, referenceFile, decodedFile, ...
%                        width, height, numFrames, blockSize, searchRange, ...
%                        I_Period, lambda, VBSEnable, FMEEnable, FastME, ...
%                        nRefFrames, QPs);

lambdatest(filename, outputFile, paddedOutputFile, referenceFile, decodedFile, ...
                       width, height, numFrames, blockSize, searchRange, ...
                       I_Period, lambda, VBSEnable, FMEEnable, FastME, ...
                       nRefFrames, QPs);

%analyze_reference_frames();

