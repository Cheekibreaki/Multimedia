addpath('../Utils');  % For utils functions
addpath('../Outputs');  % For utils functions
addpath('../Decoder');  % For Decoder functions
% Parameters
filename = '../QCIF.yuv';  % YUV file to read
% filename = '../CIF.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only_foreman.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y_foreman.yuv'; % File to store padded Y components
referenceFile = '../Outputs/referenceFrames.yuv';
decodedFile = '../Outputs/decoded_Y_foreman.yuv';
% width = 352;                     % Frame width
% height = 288;                    % Frame height
width = 176;                     % Frame width
height = 144;                    % Frame height
numFrames = 21;                 % Number of frames tfalseo process
searchRange = 16;                 % Search range r = 1,4, and 8
QP = 0;
j = 4;
VBSEnable = true;
FMEEnable = true;
FastME = true;
RCflag = 1;

if(VBSEnable == true)
    j = j-1;
end

blockSize = 2^j;                   % Block size for motion estimation
dct_blockSize = 2^j;
I_Period = 4; 
nRefFrames = 1;                 % Can take value from 1 to 4

%targetBR = 2737152;
targetBR = 1094860;


fps = 30;
%bitCountPerRow = [2112, 1520, 1256, 1156, 1076, 1056, 65, 32, 15, 7];
% For CIF
% i_bitCountPerRow = [25051, 25360, 18853, 12997, 8962, 5237, 3086, 1714, 737, 439,389,391];
% p_bitCountPerRow = [19646,14576,9974,6222,3538,2058,1293,734,475,378,354];
% For QCIF
i_bitCountPerRow = [13781,10479,7783,5642,3870,2479,1326,489,216,186,181];
p_bitCountPerRow = [9624,7059,4936,3315,2121,1318,824,433,253,189,183];

function lambda = get_lambda_for_qp(QP)
    if QP == 0
        lambda = 0.01;
    elseif QP == 1
        lambda = 0.02;
    elseif QP == 2
        lambda = 0.1667;
    elseif QP == 3
        lambda = 0.1667;
    elseif QP == 4
        lambda = 0.1667;
    elseif QP == 5
        lambda = 0.1852;
    elseif QP == 6
        lambda = 0.2037
    elseif QP == 7
        lambda = 0.222;
    elseif QP == 8
        lambda = 0.3148;
    elseif QP == 9
        lambda = 0.4074
    elseif QP == 10
        lambda = 0.5;
    elseif QP == 11
        lambda = 0.55;
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
per_block_row_budget = 0;
if RCflag == true
    per_block_row_budget = get_per_row_budget(targetBR,fps,paddedWidth,paddedHeight,blockSize);
end


% % encoder
encoder(referenceFile, paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable, FMEEnable,FastME,RCflag,per_block_row_budget,p_bitCountPerRow,i_bitCountPerRow  );
[total_byte,bytes_list] = decoder(decodedFile);
% %decoder
compareYUVFrames(referenceFile, outputFile, decodedFile, width, height, numFrames);
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

