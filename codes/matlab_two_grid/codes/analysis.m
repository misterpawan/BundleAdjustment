clc, clear
filename = 'JTJ49';
addpath(genpath('/home/pawan/work/datasets')); % from IIIT

A = load(strcat(filename, '.txt'));
A = spconvert(A);
%A = Problem.A; 
A = A + triu(A,1)';
nc = '100';
[m,n] = size(A)

P = zeros(m,1); P(1) = 1;

Ac = P'* A * P;

B = P * (Ac \ P');

