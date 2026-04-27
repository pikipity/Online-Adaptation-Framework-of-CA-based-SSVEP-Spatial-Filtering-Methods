clear;clc;
data_path = 'data/filtered_data';
save_path = 'test_online_mscca_res';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
fs = 250;
latencyDelay = round(0.14*fs);
subj_num = 35;
ch_used=[48 54 55 56 57 58 61 62 63];
ch_num = length(ch_used);
block_num = 6;
load(fullfile(data_path, 'Freq_Phase.mat'),'freqs','phases')
stim_num = length(freqs);
subband_num = 5;
fb_coef=[1:subband_num].^(-1.25)+0.25;
harmonic_num = 5;
n_neighbor = 12;
%
seed_num=42;
s = RandStream('mt19937ar','Seed',seed_num);
% p = parpool(3);
for run_ind =1:20
    %
    exp_seq = [];
    for block_ind = 1:block_num
        for stim_ind = 1:stim_num
            exp_seq = [exp_seq;
                       [stim_ind block_ind]];
        end
    end
    exp_seq = exp_seq(randperm(s,size(exp_seq, 1)),:);
    %
    for time_length = 1.5:-0.1:0.5 % 0.5:0.1:1.5 % 
        sig_len = round(time_length*fs);
        check_corr = zeros(18, subj_num, 8, size(exp_seq,1));
        check_r_diff = zeros(18, subj_num, 8, size(exp_seq,1), 3);
        for sub_ind = 1:subj_num
            load(fullfile(data_path, sprintf('S%d.mat',sub_ind)))
            % stim_num, block_num, subband_num, ch_num, sig_len
            x = eeg_data(:, :, :, :, latencyDelay+1:latencyDelay+sig_len);
            classifier = {};
            for subband_ind = 1:subband_num
                for classifier_ind = 1:18
                    classifier{classifier_ind, subband_ind} = online_mscca(freqs, phases, harmonic_num, fs,...
                                                                           ch_num, sig_len);
                end
            end
            for trial_ind = 1:size(exp_seq,1)
                test_x = squeeze(x(exp_seq(trial_ind,1), exp_seq(trial_ind,2),:,:,:));
                for classifier_ind = 1:18
                    r_1 = [];
                    r_2 = [];
                    r_3 = [];
                    r_4 = [];
                    r_5 = [];
                    r_6 = [];
                    r_7 = [];
                    r_8 = [];
                    for subband_ind = 1:subband_num
                        test_x_subband = squeeze(test_x(subband_ind, :, :)).';
                        [r1_seq, r2_seq, r3_seq, r4_seq, r5_seq, r6_seq] = classifier{classifier_ind, subband_ind}.eval(test_x_subband);
                        r_1 = [r_1;
                               r1_seq];
                        r_2 = [r_2;
                               sign(r1_seq).*r1_seq.^2 + sign(r2_seq).*r2_seq.^2];
                        r_3 = [r_3;
                               sign(r1_seq).*r1_seq.^2 + sign(r2_seq).*r2_seq.^2 + sign(r3_seq).*r3_seq.^2];
                        r_4 = [r_4;
                               sign(r1_seq).*r1_seq.^2 + sign(r2_seq).*r2_seq.^2 + sign(r3_seq).*r3_seq.^2 + sign(r4_seq).*r4_seq.^2];
                        r_5 = [r_5;
                               sign(r1_seq).*r1_seq.^2 + sign(r2_seq).*r2_seq.^2 + sign(r4_seq).*r4_seq.^2];
                        r_6 = [r_6;
                               sign(r1_seq).*r1_seq.^2 + sign(r2_seq).*r2_seq.^2 + sign(r3_seq).*r3_seq.^2 + sign(r4_seq).*r4_seq.^2 + sign(r5_seq).*r5_seq.^2 + sign(r6_seq).*r6_seq.^2];
                        if sign(r5_seq).*r5_seq.^2 + sign(r6_seq).*r6_seq.^2 == 0
                            r_7 = [r_7;
                                    r1_seq];
                        else
                            r_7 = [r_7;
                                   sign(r5_seq).*r5_seq.^2 + sign(r6_seq).*r6_seq.^2];
                        end
                        r_8 = [r_8;
                               sign(r1_seq).*r1_seq.^2 + sign(r2_seq).*r2_seq.^2 + sign(r4_seq).*r4_seq.^2 + sign(r5_seq).*r5_seq.^2 + sign(r6_seq).*r6_seq.^2];
                    end
                    r_1 = fb_coef * r_1;
                    r_2 = fb_coef * r_2;
                    r_3 = fb_coef * r_3;
                    r_4 = fb_coef * r_4;
                    r_5 = fb_coef * r_5;
                    r_6 = fb_coef * r_6;
                    r_7 = fb_coef * r_7;
                    r_8 = fb_coef * r_8;
                    for check_r_diff_ind = 1:3
                        r_sort = sort(r_1,'descend');
                        check_r_diff(classifier_ind, sub_ind, 1, trial_ind, check_r_diff_ind) = r_sort(1) - r_sort(1+check_r_diff_ind);
                        r_sort = sort(r_2,'descend');
                        check_r_diff(classifier_ind, sub_ind, 2, trial_ind, check_r_diff_ind) = r_sort(1) - r_sort(1+check_r_diff_ind);
                        r_sort = sort(r_3,'descend');
                        check_r_diff(classifier_ind, sub_ind, 3, trial_ind, check_r_diff_ind) = r_sort(1) - r_sort(1+check_r_diff_ind);
                        r_sort = sort(r_4,'descend');
                        check_r_diff(classifier_ind, sub_ind, 4, trial_ind, check_r_diff_ind) = r_sort(1) - r_sort(1+check_r_diff_ind);
                        r_sort = sort(r_5,'descend');
                        check_r_diff(classifier_ind, sub_ind, 5, trial_ind, check_r_diff_ind) = r_sort(1) - r_sort(1+check_r_diff_ind);
                        r_sort = sort(r_6,'descend');
                        check_r_diff(classifier_ind, sub_ind, 6, trial_ind, check_r_diff_ind) = r_sort(1) - r_sort(1+check_r_diff_ind);
                        r_sort = sort(r_7,'descend');
                        check_r_diff(classifier_ind, sub_ind, 7, trial_ind, check_r_diff_ind) = r_sort(1) - r_sort(1+check_r_diff_ind);
                        r_sort = sort(r_8,'descend');
                        check_r_diff(classifier_ind, sub_ind, 8, trial_ind, check_r_diff_ind) = r_sort(1) - r_sort(1+check_r_diff_ind);
                    end
                    [~,I_1] = max(r_1);
                    [~,I_2] = max(r_2);
                    [~,I_3] = max(r_3);
                    [~,I_4] = max(r_4);
                    [~,I_5] = max(r_5);
                    [~,I_6] = max(r_6);
                    [~,I_7] = max(r_7);
                    [~,I_8] = max(r_8);
                    if I_1 == exp_seq(trial_ind,1)
                        check_corr(classifier_ind, sub_ind, 1, trial_ind) = 1;
                    end
                    if I_2 == exp_seq(trial_ind,1)
                        check_corr(classifier_ind, sub_ind, 2, trial_ind) = 1;
                    end
                    if I_3 == exp_seq(trial_ind,1)
                        check_corr(classifier_ind, sub_ind, 3, trial_ind) = 1;
                    end
                    if I_4 == exp_seq(trial_ind,1)
                        check_corr(classifier_ind, sub_ind, 4, trial_ind) = 1;
                    end
                    if I_5 == exp_seq(trial_ind,1)
                        check_corr(classifier_ind, sub_ind, 5, trial_ind) = 1;
                    end
                    if I_6 == exp_seq(trial_ind,1)
                        check_corr(classifier_ind, sub_ind, 6, trial_ind) = 1;
                    end
                    if I_7 == exp_seq(trial_ind,1)
                        check_corr(classifier_ind, sub_ind, 7, trial_ind) = 1;
                    end
                    if I_8 == exp_seq(trial_ind,1)
                        check_corr(classifier_ind, sub_ind, 8, trial_ind) = 1;
                    end
                    for subband_ind = 1:subband_num
                        test_x_subband = squeeze(test_x(subband_ind, :, :)).';
                        if classifier_ind == 1 
                            classifier{classifier_ind, subband_ind}.update_parameters(test_x_subband, exp_seq(trial_ind,1))
                        elseif classifier_ind == 2
                            classifier{classifier_ind, subband_ind}.update_parameters(test_x_subband, I_1)
                        elseif classifier_ind == 3
                            classifier{classifier_ind, subband_ind}.update_parameters(test_x_subband, I_2)
                        elseif classifier_ind == 4
                            classifier{classifier_ind, subband_ind}.update_parameters(test_x_subband, I_3)
                        elseif classifier_ind == 5
                            classifier{classifier_ind, subband_ind}.update_parameters(test_x_subband, I_4)
                        elseif classifier_ind == 6
                            classifier{classifier_ind, subband_ind}.update_parameters(test_x_subband, I_5)
                        elseif classifier_ind == 7 || classifier_ind == 9 || classifier_ind == 11 || classifier_ind == 13
                            classifier{classifier_ind, subband_ind}.update_parameters(test_x_subband, I_6)
                        elseif classifier_ind == 8 || classifier_ind == 10 || classifier_ind == 12 || classifier_ind == 14
                            classifier{classifier_ind, subband_ind}.update_parameters(test_x_subband, I_7)
                        elseif classifier_ind == 15 || classifier_ind == 16 || classifier_ind == 17 || classifier_ind == 18
                            classifier{classifier_ind, subband_ind}.update_parameters(test_x_subband, I_8)
                        end
                        if classifier_ind<=8 || classifier_ind==15
                            classifier{classifier_ind, subband_ind}.update_ir('none', 0.2, [], []) % update_ir('Tikhonov', 0.2, [], []) % 
                        elseif classifier_ind==9 || classifier_ind==10 || classifier_ind==16
                            classifier{classifier_ind, subband_ind}.update_ir('Tikhonov', 0.2, [], []) % update_ir('Tikhonov', 0.2, [], []) % 
                        elseif classifier_ind==11 || classifier_ind==12 || classifier_ind==17
                            classifier{classifier_ind, subband_ind}.update_ir('ridge', 0.2, [], []) % update_ir('Tikhonov', 0.2, [], []) % 
                        elseif classifier_ind==13 || classifier_ind==14 || classifier_ind==18
                            classifier{classifier_ind, subband_ind}.update_ir('ols', 0.2, [], []) % update_ir('Tikhonov', 0.2, [], []) % 
                        end
                        classifier{classifier_ind, subband_ind}.update_gen_ssvep()
                    end
                end
                acc = [];
                for classifier_ind = 1:18
                    for class_res_type_ind = 1:8
                        acc(classifier_ind, class_res_type_ind) = mean(squeeze(check_corr(classifier_ind, sub_ind, class_res_type_ind, 1:trial_ind)));
                    end
                end
                fprintf('S%d, Time: %.1f, Trail %d, Run %d:\n', sub_ind, time_length, trial_ind, run_ind)
                disp(acc)
            end
            % save space
            for classifier_ind = 1:18
                for subband_ind = 1:subband_num
                    classifier{classifier_ind, subband_ind}.save_space();
                end
            end
            save(fullfile(save_path, sprintf('S%d_T%d_RUN%d_seed%d_classifier.mat', sub_ind, time_length*10, run_ind, seed_num)), 'classifier')
        end
        acc = [];
        for sub_ind = 1:subj_num
            for classifier_ind = 1:18
                for class_res_type_ind = 1:8
                    acc(sub_ind, classifier_ind, class_res_type_ind) = mean(squeeze(check_corr(classifier_ind, sub_ind, class_res_type_ind, :)));
                end
            end
        end
        save(fullfile(save_path, sprintf('acc_res_T%d_RUN%d_seed%d.mat', time_length*10, run_ind, seed_num)), 'acc', 'check_r_diff_ind', 'check_corr')
    end
end
% delete(p)