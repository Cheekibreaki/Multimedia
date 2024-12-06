% Clear workspace and command window
clear; clc; close all;

% --- Load Data ---
% Replace 'file1.mat' and 'file2.mat' with your actual file names.
% Replace 'A' with the actual variable name inside the .mat files if different.
data1 = load('pred1.mat'); 
data2 = load('pred2.mat'); 

% Extract the matrices (assuming they are stored under variable 'A')
A = data1.predictedFrame;  % A 288x352 double
B = data2.predictedFrame;  % B 288x352 double

% --- Check Dimensions ---
if ~isequal(size(A), size(B))
    error('The two matrices do not have the same dimensions.');
end

% --- Direct Equality Check ---
areEqual = isequal(A, B);
if areEqual
    disp('The two matrices are identical.');
else
    disp('The two matrices are not identical.');
end

% --- Compute Differences ---
D = A - B;

% Maximum absolute difference
maxDiff = max(abs(D(:)));
disp(['Maximum absolute difference: ', num2str(maxDiff)]);

% Mean difference
meanDiff = mean(D(:));
disp(['Mean difference: ', num2str(meanDiff)]);

% Root-mean-square difference
rmsDiff = sqrt(mean(D(:).^2));
disp(['RMS difference: ', num2str(rmsDiff)]);

% --- Statistical Relationship (Correlation) ---
% Flatten the matrices into vectors for correlation
corrValue = corrcoef(A(:), B(:));
% corrcoef returns a 2x2 matrix where off-diagonal elements are the correlation
disp(['Correlation between A and B: ', num2str(corrValue(1,2))]);

