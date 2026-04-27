classdef online_mstrca < handle
    properties
        ref_sig = 0;
%         Omega_Y = 0;
        Omega_X = 0;
        F = 0;
        K_cov = 0;
        W_trca = nan;
%         V_cca = nan;
        combine_matrix_F = nan;
%         combine_matrix_K = nan;
        stim_num = 0;

        com_H = 0;
        com_H_ext = 0;
        w_tlcca = nan;
        ir_diffch = 0;
        ir_filtered = 0;
        fs = 0;
        ch_num = 0;
        gen_Omega_X = 0;
        gen_template = 0;
        remove_points = [];
    end

    methods

        function obj = online_mstrca(freqs, phases, harmonic_num, fs,...
                                    ch_num, sig_len)
            L = sig_len;
            t = linspace(0, (L-1)/fs, L);
            stim_num = length(freqs);
            ref_sig = zeros(stim_num, harmonic_num*2, L);
            for stim_ind = 1:stim_num
                tmp = [];
                for harmonic_ind = 1:harmonic_num
                    tmp = [tmp;
                           sin(2*pi*harmonic_ind*freqs(stim_ind)*t+harmonic_ind*phases(stim_ind));
                           cos(2*pi*harmonic_ind*freqs(stim_ind)*t+harmonic_ind*phases(stim_ind))];
                end
                ref_sig(stim_ind,:,:) = tmp;
            end
            obj.ref_sig = permute(ref_sig, [1 3 2]);

            obj.com_H = {};
            obj.com_H_ext = {};
            for stim_ind = 1:stim_num
                obj.com_H{stim_ind} = gen_common_H(freqs(stim_ind), phases(stim_ind), fs, sig_len, 60);
                obj.com_H_ext{stim_ind} = gen_common_H(freqs(stim_ind), phases(stim_ind), fs, sig_len+floor(1/freqs(stim_ind)*fs), 60);
                obj.remove_points(stim_ind) = floor(1/freqs(stim_ind)*fs);
            end

            obj.Omega_X = zeros(stim_num, sig_len, ch_num);
            obj.F = zeros(ch_num, ch_num);
            obj.K_cov = zeros(ch_num, ch_num);
            obj.combine_matrix_F = [zeros(sig_len), eye(sig_len);
                                    eye(sig_len),   eye(sig_len)];
            obj.stim_num = stim_num;

            obj.fs = fs;
            obj.ch_num = ch_num;
            obj.gen_Omega_X = zeros(stim_num, sig_len, ch_num);
        end

        function ref_sig = get_ref_sig(obj, ind)
            ref_sig = squeeze(obj.ref_sig(ind, :, :));
        end

        function Omega_X = get_Omega_X(obj, ind)
            Omega_X = squeeze(obj.Omega_X(ind, :, :));
        end

        function Omega_X = get_gen_Omega_X(obj, ind)
            Omega_X = squeeze(obj.gen_Omega_X(ind, :, :));
        end

        function update_ir(obj, reg_method, lambda, iter_err_threshold, max_iter)
            comb_ssvep = [];
            comb_com_H = [];
            for stim_ind = 1:obj.stim_num
                tmp_X = obj.get_Omega_X(stim_ind);
                if sum(sum(tmp_X)) == 0
                else
                    comb_ssvep = [comb_ssvep tmp_X.'];
                    comb_com_H = [comb_com_H obj.com_H{stim_ind}];
                end
            end
            
            % if isnan(obj.W_cca)
            %     init_w = [];
            % else
            %     init_w = obj.W_cca(:,1).';
            % end
            init_w = [];
            [obj.ir_filtered, obj.w_tlcca, ~] = cal_ir(comb_ssvep, comb_com_H, init_w, reg_method, lambda, obj.fs, iter_err_threshold, max_iter);
            obj.w_tlcca = obj.w_tlcca.';

            obj.ir_diffch = [];
            for ch_ind = 1:obj.ch_num
                [obj.ir_diffch(ch_ind,:), ~, ~] = cal_ir(comb_ssvep(ch_ind,:), comb_com_H, [], reg_method, lambda, obj.fs, [], []);
            end
        end

        function update_gen_ssvep(obj)
            obj.gen_Omega_X = [];
            obj.gen_template = [];
            for stim_ind = 1:obj.stim_num
                for ch_ind = 1:obj.ch_num
                     tmp = obj.ir_diffch(ch_ind,:) * obj.com_H_ext{stim_ind};
                     obj.gen_Omega_X(stim_ind,:,ch_ind) = tmp(obj.remove_points(stim_ind)+1:end);
                end
                 tmp = obj.ir_filtered * obj.com_H_ext{stim_ind};
                 obj.gen_template(stim_ind,:) = tmp(obj.remove_points(stim_ind)+1:end);
            end
        end

        function update_parameters(obj, X_new, predict_ind)
            tmp = [obj.get_Omega_X(predict_ind); % Omega_X (k)
                   X_new];
            obj.F = obj.F + tmp.' * obj.combine_matrix_F * tmp; % F (k+1)

            obj.K_cov = obj.K_cov + X_new.'*X_new; % K_cov (k+1)

            obj.Omega_X(predict_ind, :, :) = obj.get_Omega_X(predict_ind) + X_new; % Omega_X (k+1)

            [W_store, ~] = solve_lsframe(obj.F, obj.K_cov);
            obj.W_trca = real(W_store{end});
        end

        function [r1_seq, r3_seq, r4_seq, r5_seq, r6_seq] = eval(obj, X)
            r1_seq = zeros(1, obj.stim_num);
%             r2_seq = r1_seq;
            r3_seq = r1_seq;
            r4_seq = r1_seq;
            r5_seq = r1_seq;
            r6_seq = r1_seq;

            for stim_ind = 1:obj.stim_num
                [~,~,r_tmp] = canoncorr(X, obj.get_ref_sig(stim_ind));
                r1_seq(stim_ind) = r_tmp(1);
            end
            if isnan(obj.W_trca)
            else
                for stim_ind = 1:obj.stim_num
%                     r_tmp = corrcoef(X * obj.W_cca(:,1), obj.get_ref_sig(stim_ind) * obj.V_cca(:,1));
%                     r2_seq(stim_ind) = r_tmp(1,2);
                    r_tmp = corrcoef(X * obj.W_trca(:,1), (obj.get_gen_Omega_X(stim_ind) + obj.get_Omega_X(stim_ind)) * obj.W_trca(:,1));
                    if isnan(r_tmp(1,2))
                        r3_seq = zeros(1, obj.stim_num);
                    else
                        r3_seq(stim_ind) = r_tmp(1,2);
                    end
                    [~,~,r_tmp] = canoncorr(X * obj.W_trca(:,1), obj.get_ref_sig(stim_ind));
                    r4_seq(stim_ind) = r_tmp(1);
                end
            end
            if isnan(obj.w_tlcca)
            else
                for stim_ind = 1:obj.stim_num
                    r_tmp = corrcoef(X * obj.w_tlcca(:,1), obj.gen_template(stim_ind,:).');
                    r5_seq(stim_ind) = r_tmp(1,2);
                    [~,~,r_tmp] = canoncorr(X * obj.w_tlcca(:,1), obj.get_ref_sig(stim_ind));
                    r6_seq(stim_ind) = r_tmp(1);
                end
            end
            
        end

        function save_space(obj)
            obj.ref_sig = nan;
        end


    end
end