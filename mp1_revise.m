%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECE 417 - MP 1
% Yucheng Liang     2018@UIUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparation
clc,clear,close all
[y,Fs] = audioread('s5.wav');
[leng,~] = size(y);
frameLength = 25; % in ms
frameShift = 10; % in ms
K = frameShift/1000*Fs; % frame 'skip'
L = frameLength /1000*Fs; % frame length
fnum = leng/K-2;
clear frameLength frameShift
raw_data = zeros(L,fnum);
for count = 1:fnum
    start = (count-1)*K+1;
    raw_data(:,count) = y(start:start+L-1);
end
clear count leng start

%% spectrogram
spectrogram(y,hamming(L),(L-K),1024);

%% Analysis - Pitch Estimation
% Find P0 %%%%%%%%%%%%%%%%%%%%
% hamming window
w = hamming(L);
w_data = repmat(w,[1,fnum]).*raw_data; % windowed data
% fourier transform
kw_data = fft(w_data); %Fn_x2k(w_data,1);
% autocorrelation
phi = zeros(fnum, 2*L-1);
P0 = zeros(fnum,1);
parfor num = 1:(fnum)
    temp = w_data(:,num).*hamming(L);
    phi(num,:) = transpose(conv(temp,flip(temp,1)));
    P0(num,:) = opt(phi(num,:));
end
clear num temp

%% analysis - Pitch Refinement and Spectral Envelope
% Find Am %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Am = zeros(fnum,92);
P = zeros(size(P0));
dec = Am;
tic
for n = 1:fnum
    error = Inf;
    original_pn = P0(n);
    E_v = zeros(1,1024);
    for P_n = (original_pn-2):0.2:(original_pn+2)
        esum = 0;
        Atemp = zeros(1,92);
        dtemp = zeros(1,92);
        for m = 1:(P_n-1)
            % preparation
            Sw = fft(raw_data(:,n)); %Fn_x2k(raw_data(:,n),1);
            am = ceil((m-1/2)/P_n*L);
            bm = floor((m+1/2)/P_n*L);
            % voiced
            shift_num = ceil(L/2-m/P_n*L);
            E_v = fft(circshift(hamming(L),-shift_num)); %Fn_x2k(hamming(L),1); % shift
            temp1 = Sw.*conj(E_v);
            temp1 = sum(temp1(am:bm));
            temp2 = abs(E_v).*abs(E_v);
            temp2 = sum(temp2(am:bm));
            A_v = temp1/temp2;
            e_v = abs(Sw-A_v.*E_v).*abs(Sw-A_v.*E_v);
            e_v = 1/2/pi*sum(e_v(am:bm));
            % unvoiced
            E_u = ones(200,1);
            temp1 = Sw.*conj(E_u);
            temp1 = sum(temp1(am:bm));
            temp2 = abs(E_u).*abs(E_u);
            temp2 = sum(temp2(am:bm));
            A_u = temp1/temp2;
            e_u = abs(Sw-A_u.*E_u).*abs(Sw-A_u.*E_u);
            e_u = 1/2/pi*sum(e_u(am:bm));
            % compare error
            if (e_v < e_u)
                Atemp(m) = A_v;
                em = e_v;
                dtemp(m) = 2; % voiced for 2
            else
                Atemp(m) = A_u;
                em = e_u;
                dtemp(m) = 1; % unvoiced for 1
            end
            esum = esum + em;
        end
        if (esum<error)
            error = esum;
            Am(n,:) = Atemp;
            P(n) = P_n;
            dec(n,:) = dtemp;
        end
    end
end
toc
clear w w_data

%% 
Am = real(Am);
s = zeros(size(y));
for f = 1:(fnum-1)
    dec_f = dec(f,:);
    Pf = P(f);
    Pf1 = P(f+1,1);
    n = transpose(1:K);
    n = n+f*K;
    for band = 1:(floor(Pf)-1)
        if (dec_f(band)==1) %1 for unvoiced
            am = ceil((band-1/2)/Pf*L);
            bm = floor((band+1/2)/Pf*L);
            Sw = Fn_x2k(raw_data(:,band),1);
            var = 1/(bm-am)*sum(sqrt(Sw));
            U = randn(K,1)*var + 1i*randn(K,1)*var;
            uf = Fn_k2x(U,1);
            s(f*K+1:(f+1)*K) = s(f*K+1:(f+1)*K) + (f+1-n/K).*uf;
            if f~=1
                s((f-1)*K+1:f*K) = s((f-1)*K+1:f*K) + (n/K-f).*uf;
            end
        elseif (dec_f(band)==2) % 2 for voiced
            w0 = (f+1-n/K)*2*pi/Pf+(n/K-f)*2*pi/Pf1;
            theta_m = zeros(80,1);
            for ind = 2:80   % get theta
                theta_m(ind) = theta_m(ind-1)+band*w0(ind);
            end
            if (dec_f(band+1)~=2)   % check following
                A_m = (f+1-n/K)*Am(f,band);
                s(f*K+1:(f+1)*K) = s(f*K+1:(f+1)*K) + A_m.*cos(theta_m);
            end
            if (band~=1)    % check first
                n_temp = n-K;
                A_mf = 0;
                if (dec_f(band-1)==2)
                    A_mf = Am(f,band-1);
                end
                A_m = (f+1-n_temp/K)*A_mf+(n_temp/K-f)*Am(f,band);
                s((f-1)*K+1:f*K) = s((f-1)*K+1:f*K) + A_m.*cos(theta_m);
            end
        end   
    end
end

%%
audiowrite('trail1.wav',real(s./(max(abs(s)))),Fs);

%% step3

for f = 1:(fnum-1)
    dec_f = dec(f,:);
    for 
    if (dec_f(



