clear;clc;
data_path = 'data';
save_path = 'data/filtered_data';
if ~exist(save_path)
    mkdir(save_path)
end
Fs = 250;
% latencyDelay = round(0.14*Fs);
subj_num = 35;
ch_used=[48 54 55 56 57 58 61 62 63];
ch_num = length(ch_used);
block_num = 6;
load(fullfile(data_path, 'Freq_Phase.mat'),'freqs','phases')
stim_num = length(freqs);
subband_num = 5;
sig_len = round(5*Fs);
% Build filter bank
subband_signal = [];
for k = 1:subband_num
    Wp = [(8*k)/(Fs/2) 90/(Fs/2)];
    Ws = [(8*k-2)/(Fs/2) 100/(Fs/2)];
    [N,Wn] = cheb1ord(Wp,Ws,3,40);
    [subband_signal(k).bpB,subband_signal(k).bpA] = cheby1(N,0.5,Wn);
end
% Build notch filter
Fo = 50;
Q = 35;
BW = (Fo/(Fs/2))/Q;
[notchB,notchA] = iircomb(Fs/Fo,BW,'notch');
% Filter signal
for sub_ind = 1:subj_num
    fprintf('Filter S%d\n',sub_ind)
    load(fullfile(data_path,sprintf('S%d.mat',sub_ind)))
    data = permute(data,[3,4,1,2]);
    data = data(:,:,ch_used,(floor(0.5*Fs)+1):(floor(0.5*Fs)+sig_len));
    eeg_data = zeros(stim_num, block_num, subband_num, ch_num, sig_len);
    for stim_ind = 1:stim_num
        for block_ind = 1:block_num
            for subband_ind = 1:subband_num
                y0 = squeeze(data(stim_ind, block_ind, :, :));
                y1 = filtfilt(notchB, notchA, y0.').';
                y2 = filtfilt(subband_signal(subband_ind).bpB,subband_signal(subband_ind).bpA,y1.').';
                eeg_data(stim_ind, block_ind, subband_ind, :, :) = y2;
            end
        end
    end
    save(fullfile(save_path, sprintf('S%d.mat',sub_ind)), 'eeg_data')
end