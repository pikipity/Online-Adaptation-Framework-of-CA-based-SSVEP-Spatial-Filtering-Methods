function com_H = gen_common_H(f0, p, Fs, sig_len, varargin)
% com_H = gen_common_H(f0, p, Fs, tw, refresh_rate)
% refresh_rate: default is 60
    
    if nargin<4
        error('Input variables are not enough.')
    elseif nargin == 4
        refresh_rate = 60;
    elseif nargin == 5
        refresh_rate = varargin{1};
    else
        error('Input variables are too many.')
    end
    tw = sig_len/Fs;
    [H, ~, ~] = conv_H_0725(f0,p,Fs,tw,refresh_rate);
    [com_H] = common_H(H, f0);

end