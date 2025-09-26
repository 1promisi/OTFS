clc
clear all
close all
tic
%% OTFS parameters%%%%%%%%%%
% number of symbol
N = 32;
% number of subcarriers
M = 64;
% size of constellation
M_mod = 2;
M_bits = log2(M_mod);
B = 4*10^3;
delta_f = B / M;       % 子载波间隔
T = 1 / delta_f;          % OTFS符号周期
sampling_rate = M * delta_f;
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
% number of symbols per frame
N_syms_perfram = N*M;
% number of bits per frame 
N_bits_perfram = N*M*M_bits;

SNR_dB = 0:5:20;
SNR = 10.^(SNR_dB/10);
noise_var_sqrt = sqrt(1./SNR);
sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;

rng(1)
N_fram = 5;
err_ber = zeros(length(SNR_dB),1);

for iesn0 = 1:length(SNR_dB)
    for ifram = 1:N_fram
        %% random input bits generation%%%%%
        data_info_bit = randi([0,1],N_bits_perfram,1);
        data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
        x = qammod(data_temp,M_mod,'gray');
        x = reshape(x,N,M);
        
        %% OTFS modulation%%%%
        % ISFFT
        X = fft(ifft(x).').'/sqrt(M/N); 
        % Heisenberg transform
        s_mat = ifft(X.')*sqrt(M); 
        s = s_mat(:);
        
        %% OTFS channel generation%%%%
%         taps = 3;
%         delay_taps = [0 1 4];
%         Doppler_taps = [0 -1 2];
%         pdp = [0 -1.0 -2.0];
%         pow_prof = 10.^(pdp/10);
%         pow_prof = pow_prof/sum(pow_prof);%normalization of power delay profile
%         chan_coef = sqrt(pow_prof).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));%channel coef. for each path
        tmp = load('C:\Users\86198\Desktop\mobile receiver\OTFS1\sample_OTFS\otfs_channel_struct_delta_tau.mat'); 
        data_cell = tmp.impulse_responses;

        num_steps = 1;
        max_taps = 0;
        for n = 1:num_steps
            max_taps = max(max_taps, data_cell{n}.num_echoes);
        end
        taps = max_taps;

        delay_taps_all   = zeros(taps, num_steps);
        Doppler_taps_all = zeros(taps, num_steps);
        chan_coef_all    = zeros(taps, num_steps);

        pdp = [0 -1.0 -2.0 -3.0 -8.0 -17.2];  
        if length(pdp) < taps
            pdp = [pdp, repmat(pdp(end), 1, taps - length(pdp))];
        end
        pow_prof = 10.^(pdp/10);
        pow_prof = pow_prof / sum(pow_prof);

        for n = 1:num_steps
            s_ch = data_cell{n}; 
            num_echoes = s_ch.num_echoes;
            delay_s   = zeros(taps, 1);
            phase_s   = zeros(taps, 1);
            
            if num_echoes > 0
                delay_s(1:num_echoes) = s_ch.delay(1:min(num_echoes, end));
                phase_s(1:num_echoes) = s_ch.Doppler_shift(1:min(num_echoes, end));
            end
            
            delay_taps = round(delay_s * 1e3 / sampling_rate);  
            Doppler_taps = round(phase_s / (N * T));
            
            delay_taps_all(:, n)   = delay_taps;
            Doppler_taps_all(:, n) = Doppler_taps;
            
            pow_weights = sqrt(pow_prof(1:taps))';
            if size(pow_weights, 2) > 1
                pow_weights = pow_weights';
            end
           
            chan_coef = zeros(taps, 1);
            if num_echoes > 0
                random_complex = (randn(num_echoes, 1) + 1i*randn(num_echoes, 1));
                chan_coef(1:num_echoes) = pow_weights(1:num_echoes) .* (sqrt(1/2) * random_complex);
            end
            chan_coef_all(:, n) = chan_coef;
        end

        delay_taps = delay_taps_all;
        Doppler_taps = Doppler_taps_all;
        chan_coef = chan_coef_all;

        
        %% OTFS channel output%%%%%
        L = max(delay_taps);
        s_cp = [s(N*M-L+1:N*M);s]; % add CP
        s_chan = 0;
        for itao = 1:taps
            s_chan = s_chan+chan_coef(itao)*circshift([s_cp.*exp(1j*2*pi/M ...
                *(-L:-L+length(s_cp)-1)*Doppler_taps(itao)/N).';zeros(delay_taps(end),1)],delay_taps(itao));
        end
        noise = sqrt(sigma_2(iesn0)/2)*(randn(size(s_chan)) + 1i*randn(size(s_chan)));
        r = s_chan + noise;
        r = r(L+1:L+(N*M)); % discard CP
        
        %% OTFS demodulation%%%%
        r_mat = reshape(r,M,N);
        Y = fft(r_mat)/sqrt(M); % Wigner transform
        Y = Y.';
        y = ifft(fft(Y).').'/sqrt(N/M); % SFFT
        
        %% message passing detector%%%%
        yv = reshape(y,N*M,1);
        n_ite = 200;
        delta_fra = 0.6;
        alphabet = qammod(0:M_mod-1,M_mod,'gray');
        mean_int = zeros(N*M,taps);
        var_int = zeros(N*M,taps);
        p_map = ones(N*M,taps,M_mod)*(1/M_mod);
        conv_rate_prev = -0.1;

        for ite=1:n_ite
            for ele1=1:1:M
                for ele2=1:1:N
                    mean_int_hat = zeros(taps,1);
                    var_int_hat = zeros(taps,1);
                    for tap_no=1:taps
                        m = ele1-1-delay_taps(tap_no)+1;
                        add_term = exp(1i*2*(pi/M)*(m-1)*(Doppler_taps(tap_no)/N));
                        add_term1 = 1;
                        if ele1-1<delay_taps(tap_no)
                            n = mod(ele2-1-Doppler_taps(tap_no),N) + 1;
                            add_term1 = exp(-1i*2*pi*((n-1)/N));
                        end
                        new_chan = add_term * (add_term1) * chan_coef(tap_no);
                        
                        for i2=1:1:M_mod
                            mean_int_hat(tap_no) = mean_int_hat(tap_no) + p_map(N*(ele1-1)+ele2,tap_no,i2) * alphabet(i2);
                            var_int_hat(tap_no) = var_int_hat(tap_no) + p_map(N*(ele1-1)+ele2,tap_no,i2) * abs(alphabet(i2))^2;
                        end
                        mean_int_hat(tap_no) = mean_int_hat(tap_no) * new_chan;
                        var_int_hat(tap_no) = var_int_hat(tap_no) * abs(new_chan)^2;
                        var_int_hat(tap_no) = var_int_hat(tap_no) - abs(mean_int_hat(tap_no))^2;
                    end
                    
                    mean_int_sum = sum(mean_int_hat);
                    var_int_sum = sum(var_int_hat)+(sigma_2(iesn0));
                    
                    for tap_no=1:taps
                        mean_int(N*(ele1-1)+ele2,tap_no) = mean_int_sum - mean_int_hat(tap_no);
                        var_int(N*(ele1-1)+ele2,tap_no) = var_int_sum - var_int_hat(tap_no);
                    end
                end
            end
            
            sum_prob_comp = zeros(N*M,M_mod);
            dum_eff_ele1 = zeros(taps,1);
            dum_eff_ele2 = zeros(taps,1);
            for ele1=1:1:M
                for ele2=1:1:N
                    dum_sum_prob = zeros(M_mod,1);
                    log_te_var = zeros(taps,M_mod);
                    for tap_no=1:taps
                        if ele1+delay_taps(tap_no)<=M
                            eff_ele1 = ele1 + delay_taps(tap_no);
                            add_term = exp(1i*2*(pi/M)*(ele1-1)*(Doppler_taps(tap_no)/N));
                            int_flag = 0;
                        else
                            eff_ele1 = ele1 + delay_taps(tap_no)- M;
                            add_term = exp(1i*2*(pi/M)*(ele1-1-M)*(Doppler_taps(tap_no)/N));
                            int_flag = 1;
                        end
                        add_term1 = 1;
                        if int_flag==1
                            add_term1 = exp(-1i*2*pi*((ele2-1)/N));
                        end
                        eff_ele2 = mod(ele2-1+Doppler_taps(tap_no),N) + 1;
                        new_chan = add_term * add_term1 * chan_coef(tap_no);
                        
                        dum_eff_ele1(tap_no) = eff_ele1;
                        dum_eff_ele2(tap_no) = eff_ele2;
                        for i2=1:1:M_mod
                            dum_sum_prob(i2) = abs(yv(N*(eff_ele1-1)+eff_ele2)- mean_int(N*(eff_ele1-1)+eff_ele2,tap_no) - new_chan * alphabet(i2))^2;
                            dum_sum_prob(i2)= -(dum_sum_prob(i2)/var_int(N*(eff_ele1-1)+eff_ele2,tap_no));
                        end
                        dum_sum = dum_sum_prob - max(dum_sum_prob);
                        dum1 = sum(exp(dum_sum));
                        log_te_var(tap_no,:) = dum_sum - log(dum1);
                    end
                    for i2=1:1:M_mod
                        ln_qi(i2) = sum(log_te_var(:,i2));
                    end
                    dum_sum = exp(ln_qi - max(ln_qi));
                    dum1 = sum(dum_sum);
                    sum_prob_comp(N*(ele1-1)+ele2,:) = dum_sum/dum1;
                    for tap_no=1:1:taps
                        eff_ele1 = dum_eff_ele1(tap_no);
                        eff_ele2 = dum_eff_ele2(tap_no);
                        dum_sum = log_te_var(tap_no,:);
                        ln_qi_loc = ln_qi - dum_sum;
                        dum_sum = exp(ln_qi_loc - max(ln_qi_loc));
                        dum1 = sum(dum_sum);
                        p_map(N*(eff_ele1-1)+eff_ele2,tap_no,:) = (dum_sum/dum1)*delta_fra + (1-delta_fra)*reshape(p_map(N*(eff_ele1-1)+eff_ele2,tap_no,:),1,M_mod);
                    end
                end
            end
            conv_rate =  sum(max(sum_prob_comp,[],2)>0.99)/(N*M);
            if conv_rate==1
                sum_prob_fin = sum_prob_comp;
                break;
            elseif conv_rate > conv_rate_prev
                conv_rate_prev = conv_rate;
                sum_prob_fin = sum_prob_comp;
            elseif (conv_rate < conv_rate_prev - 0.2) && conv_rate_prev > 0.95
                break;
            end
        end
        x_est = zeros(N,M);
        for ele1=1:1:M
            for ele2=1:1:N
                [~,pos] = max(sum_prob_fin(N*(ele1-1)+ele2,:));
                x_est(ele2,ele1) = alphabet(pos);
            end
        end

        %% output bits and errors count%%%%%
        data_demapping = qamdemod(x_est,M_mod,'gray');
        data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
        errors = sum(xor(data_info_est,data_info_bit));
        err_ber(iesn0) = errors + err_ber(iesn0);
    end
end

err_ber_fram = err_ber/N_bits_perfram./N_fram;
semilogy(SNR_dB, err_ber_fram,'-*','LineWidth',2);
title(sprintf('OTFS'))
ylabel('BER'); xlabel('SNR in dB');grid on
toc
