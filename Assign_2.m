clc
clear all
format long
HW_Data = readmatrix('Hotwire.dat');
N = size(HW_Data, 1); %Total No. of Frames
T = round(max(HW_Data(:, 1))); %Total Time Period
f = N/T; %Sampling Frequency
[autocorr, lag] = xcorr(HW_Data(:, 2) - mean(HW_Data(:,2)), 'coeff');
lag_ms = lag/f *1000; %Lag in milliseconds 
plot(lag_ms, autocorr)
title('Normalized Autocorrelation for Hot Wire Data')
xlabel('lags (in ms)')
ylabel('Correlation Coefficient')
TS = trapz(lag_ms, autocorr); 
disp(TS);
%%
