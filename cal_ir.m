function [ir, w, iter_err] = cal_ir(ssvep, com_H, init_w, reg_method, lambda, Fs, iter_err_threshold, max_iter)
    n_ch = size(ssvep,1);
    if strcmp(reg_method, 'none')
        M_H0 = 0;
        M_ssvep = 0;
    else
        M_H0 = regmat(size(com_H,1),reg_method) * lambda * Fs;
        M_ssvep = regmat(size(ssvep,1),reg_method) * lambda * Fs;
    end
    if isempty(init_w)
        init_w = randn(1,n_ch);
    end
    if isempty(iter_err_threshold)
        iter_err_threshold = 0.0001;
    end
    if isempty(max_iter)
        max_iter = 200;
    end
    
    if n_ch == 1 % 1 channel eeg does not need iteration
        ir = ssvep * com_H' * pinv(com_H * com_H' + M_H0);
        w = 1;
        iter_err = 0;
    else
        w0_old=init_w;
        x_hat_old=w0_old*ssvep*com_H'*inv(com_H*com_H' + M_H0);
        e_old=norm(w0_old*ssvep-x_hat_old*com_H);
        iter_err=100;
        iter=1;
        while (iter_err(iter)>0.0001 && iter<200)
            w0_new=x_hat_old*com_H*ssvep'*inv(ssvep*ssvep' + M_ssvep);
            x_hat_new=w0_new*ssvep*com_H'*inv(com_H*com_H' + M_H0);
            e_new=norm(w0_new*ssvep-x_hat_new*com_H);
            iter=iter+1;
            iter_err(iter)=abs(e_old-e_new);
            w0_old=w0_new;
            w0_old=w0_old/std(w0_old);
            x_hat_old=x_hat_new;
            x_hat_old=x_hat_old/std(x_hat_old);
            e_old=e_new;
        end
        ir=x_hat_new;
        w=w0_new;
    end
end