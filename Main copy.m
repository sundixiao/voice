%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECE 417 - MP 1
% Yucheng Liang     2018@UIUC
% Last change: synthesis + window shift number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparation
clc,clear,close all
tic
[y0,Fs] = audioread('s5.wav');
y = [y0;zeros(120,1)];


init_mean = sum(y)/size(y,1);
init_energy = sum(power(y,2));
y = (y-init_mean)/init_energy; 
frameLength = 25; % in ms
frameShift = 10; % in ms
K = frameShift/1000*Fs; % frame 'skip'
L = frameLength /1000*Fs; % frame length
fnum = 300;

%% Add Gaussian Noise
% var = 0.05*0.05;
% y = y + var*randn(size(y));


%% Plot--0 Spectrogram
%figure 
%spectrogram(y,hamming(L),(L-K),1024);
figure
pspectrum(y,Fs,'spectrogram')
title('Spectrogram')
%%
clear frameLength frameShift
raw_data = zeros(L,fnum);
parfor count = 1:fnum
    start = (count-1)*K+1;
    raw_data(:,count) = y(start:start+L-1);
end
clear count leng start

%% Analysis - Pitch Estimation
% Find P0 %%%%%%%%%%%%%%%%%%%%
% hamming window
w = hamming(L);
w_data = repmat(w,[1,fnum]).*raw_data; % windowed data
% fourier transform
padding = 1024;
kw_data = fft(w_data,padding); 
% autocorrelation
phi = zeros(fnum, 2*L-1);
P0 = zeros(fnum,1);
parfor num = 1:(fnum)
    temp = w_data(:,num).*hamming(L);
    phi(num,:) = transpose(conv(temp,flip(temp,1)));
    P0(num,:) = FindP(phi(num,:));
end
clear num phi temp

%% analysis - Pitch Refinement and Spectral Envelope
% Find Am %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Am = zeros(fnum,92);
P = zeros(size(P0));
dec = Am;
error = zeros(fnum);

parfor n = 1:fnum
    error(n) = Inf;
    original_pn = P0(n);
    for P_n = (original_pn-2):0.2:(original_pn+2)
        esum = 0;
        Atemp = zeros(1,92);
        dtemp = zeros(1,92);
%         ct_1 = 0;
%         ct_2 = 0;
        for m = 1:(P_n-1)
            % preparation
            Sw = fft(raw_data(:,n),padding); 
            am = ceil((m-1/2)/P_n*padding);
            bm = floor((m+1/2)/P_n*padding);
            % voiced
            shift_num = ceil(m/P_n*padding)
            E_v = circshift(fft(hamming(L),1024),shift_num);
            temp1 = Sw.*conj(E_v);
            temp1 = sum(temp1(am:bm));
            temp2 = abs(E_v).*abs(E_v);
            temp2 = sum(temp2(am:bm));
            A_v = temp1/temp2;
            err_v = abs(Sw-A_v.*E_v).*abs(Sw-A_v.*E_v);
            err_v = 1/2/pi*sum(err_v(am:bm));
            % unvoiced
            E_u = ones(padding,1);
            temp1 = Sw.*conj(E_u);
            temp1 = sum(temp1(am:bm));
            temp2 = abs(E_u).*abs(E_u);
            temp2 = sum(temp2(am:bm));
            A_u = temp1/temp2;
            err_u = abs(Sw-A_u.*E_u).*abs(Sw-A_u.*E_u);
            err_u = 1/2/pi*sum(err_u(am:bm));
            % compare error
            if (err_v <= err_u)
                Atemp(m) = A_v;
                errm = err_v;
                dtemp(m) = 2; % voiced for 2
%                 ct_2 = ct_2 + 1;
            else
                Atemp(m) = A_u;
                errm = err_u;
                dtemp(m) = 1; % unvoiced for 1
%                 ct_1 = ct_1 + 1;
            end
            esum = esum + errm;
        end
        if (esum < error(n))
            error(n) = esum;
            Am(n,:) = Atemp;
            P(n) = P_n;
%             if(ct_2<ct_1)
%                 dtemp = ones(1,92);
%             else
%                 dtemp = ones(1,92)*2;
%             end
            dec(n,:) = dtemp;
        end
    end
end
clear w w_data

%% Plot--1 Final Pitch Estimate
figure
plot(P)
title('Refined Pitch Estimate')
xlabel('Frame Number')
ylabel('Pitch Estimate')

%% Plot--2 Error per Frame
figure
plot(error)
title('Error per Frame')
xlabel('Frame Number')
ylabel('Error')





%% Synthesis
Am = real(Am);
s = zeros(K,fnum);
parfor f = 1:(fnum-1)
    dec_f = dec(f,:);
    Pf = P(f);
    Pf1 = P(f+1);
    n = transpose(1:K);
    n = n+f*K;
    s_temp = zeros(K,1);
    for band = 1:floor(Pf)-1
        % voiced
        w0 = (f+1-n/K)*2*pi/Pf+(n/K-f)*2*pi/Pf1;
        theta_m = zeros(K,1);
        theta_m(1) = band*w0(1);
        for ind = 2:K   % get theta
            theta_m(ind) = theta_m(ind-1)+band*w0(ind);
        end
        if (dec(f,band)==1)  % check Am,f and Am,f+1
            Amf = 0;
        else 
            Amf = Am(f,band);
        end
        if (dec(f+1,band)==1)
            Amf1 = 0;
        else 
            Amf1 = Am(f+1,band);
        end
        A_m = (f+1-n/K)*Amf + (n/K-f)*Amf1;
        s_temp = s_temp + A_m.*cos(theta_m);
        if (dec_f(band)==1) %1 for unvoiced
            am = ceil((band-1/2)/Pf*L);
            bm = floor((band+1/2)/Pf*L);
            Sw = fft(raw_data(:,f)); 
            var = 1/(bm-am)*sum(Sw(am:bm).*conj(Sw(am:bm)));
            U = randn(K,1)*0.5*var + 1i*randn(K,1)*0.5*var;
            uf = ifft(U); 
            s_temp = s_temp + (f+1-n/K).*uf;
            if (f>(floor(Pf)-2) || band >= floor(Pf1))
                uf1 = 0;
            else
                am1 = ceil((band-1/2)/Pf1*L);
                bm1 = floor((band+1/2)/Pf1*L);
                Sw1 = fft(raw_data(:,f+1)); 
                var1 = 1/(bm1-am1)*sum(Sw1(am1:bm1).*conj(Sw(am1:bm1)));
                U1 = randn(K,1)*0.5*var1 + 1i*randn(K,1)*0.5*var1;
                uf1 = ifft(U1); 
            end
            s_temp = s_temp + (n/K-f).*uf1;
        end  
    end
    s(:,f) = s_temp; 
end
s = reshape(s,size(y0));
s_out = real(s.*init_energy + init_mean);
s_out = s_out./(max(abs(s_out)));
toc

%%
sound(s_out,Fs);
audiowrite('trail4.wav',s_out,Fs);


figure

pspectrum(s_out,Fs,'spectrogram')
title('Spectrogram')
