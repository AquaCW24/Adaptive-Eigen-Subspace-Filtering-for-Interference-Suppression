%% Import Raw Data
% Filename: main.m
% Description: Eliminate the effects of RF interference to preserve the original signal
% Build Date: 2025/01/15
% Last edited: 2025/08/01
% Author: CHE-WEN, KUO
% Version: 1.0.0

% Revision:
% 2025/08/01: Release the official version

close all;clear;clc
%% Load path
addpath('function/')
path_folder = 'data/';
path_data = [path_folder,'RawData/'];

filenamelist = dir(path_data);
filenamelist(1:2) = [];
filelength = length(filenamelist);
%% Process

z = 1; % file no.

Header;

indd = 1:nrg1 - codelength;
data2 = reshape(data,Nchannel,fSampleLength/Nchannel);
for nc = 1:Nchannel/2
    data3 = reshape(data2(2*nc-1,:)+1j*data2(2*nc,:),nrg1,freq_dwell/pnum/cint,length(freq));
    data3 = data3(indd(1:end-2),:,1:end-1);
    data3(1,:,:) = [];
    data4 = squeeze(sum(data3,2));
end
data_rawabs = abs(data4);
data_power = data_rawabs.^2;

% Find RFI index
[CFAR_RFI_index, threshold] = CFAR(data_power, 4, 2, 0.16, "CA", "false");

[DATA_EVD, DATA_EVD_3D] = EVD(data3, CFAR_RFI_index);