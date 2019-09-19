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
s = spectrogram(y,hamming(L),(L-K),1024);

%% Analysis - Pitch Estimation
% Find P0 %%%%%%%%%%%%%%%%%%%%
% hamming window
w = hamming(L);
w_data = repmat(w,[1,fnum]).*raw_data; % windowed data
% fourier transform
kw_data = Fn_x2k(w_data,1);
% autocorrelation
phi = zeros(fnum, 2*L-1);
P0 = zeros(fnum,1);
for num = 1:(fnum)
    temp = w_data(:,num).*hamming(L);
    phi(num,:) = transpose(conv(temp,fliplr(temp)));
    P0(num,:) = opt(phi(num,:));
end
clear num phi temp

%% Pitch Refinement and Spectral Envelope
% Find Am %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Am = zeros(fnum,92);
erro_mat = zeros(fnum,92);
P_new = zeros(size(P));
tic
%errormat = zeros(P_n-1,fnum);
for n = 1:fnum
    error = 0;
    original_pn = P(n);
    for P_n = (original_pn-2):0.2:(original_pn+2)

        esum = 0;
        Atemp = zeros(1,92);
        for m = 1:(P_n-1)
            
            % preparation
            %Sw = Fn_x2k(raw_data(:,n),1);
            Sw = fft(frames_hamming(n,:));
            am = ceil((m-1/2)/P_n*L);
            bm = floor((m+1/2)/P_n*L);
            % voiced
            %E_v = Fn_x2k(hamming(L),1);
            E_v = fft(hamming(L));
            E_v = transpose(E_v);
            temp1 = Sw.*conj(E_v);
            temp1 = sum(temp1(am:bm));
            temp2 = E_v.*E_v;
            temp2 = sum(temp2(am:bm));
            A_v = temp1/temp2;
            e_v = (Sw-A_v.*E_v).*(Sw-A_v.*E_v);
            e_v = 1/2/pi*sum(e_v(am:bm));
            % unvoiced
            E_u = ones(1,200);
            temp1 = Sw.*conj(E_u);
            temp1 = sum(temp1(am:bm));
            temp2 = E_u.*E_u;
            temp2 = sum(temp2(am:bm));
            A_u = temp1/temp2;
            e_u = (Sw-A_u.*E_u).*(Sw-A_u.*E_u);
            e_u = 1/2/pi*sum(e_u(am:bm));
            % compare error
            if (e_v < e_u)
                Atemp(m) = A_v;
                em = e_v;
                erro_mat(n,m) = 2;
            else
                Atemp(m) = A_u;
                em = e_u;
                erro_mat(n,m) = 1;
            end
            esum = esum + em;
        end
        
        if (esum<error)
            error = esum;
            Am(n,:) = Atemp;
            P_new(n) = P_n;
        end
    end
end
toc
clear P0 w w_data

%%
%Unvoiced
s_u_n = zeros(1,fnum*200);
for n = 1:fnum
          for m = 1:(P(n,1)-1)
                if(erro_mat(n,m) == 1)
                am = ceil((m-1/2)/P(n,1)*L);
                bm = floor((m+1/2)/P(n,1)*L);
                Sw = Fn_x2k(raw_data(:,n),1);
                temp = Sw.*Sw;
                sum_temp = sum(abs(temp(am:bm)));
                var_sq = 1/(bm-am)*sum_temp;
                U_k = sqrt(var_sq).*randn(1,200)+1i*sqrt(var_sq).*randn(1,200);
                u_n = ifft(U_k);
                u_n = ifftshift(abs(u_n));
                end
                
          end
          for m = 1:(P(n+1,1)-1)
              if(erro_mat(n+1,m) == 1)
                am = ceil((m-1/2)/P(n+1,1)*L);
                bm = floor((m+1/2)/P(n+1,1)*L);
                Sw = Fn_x2k(raw_data(:,n+1),1);
                temp = Sw.*Sw;
                sum_temp = sum(abs(temp(am:bm)));
                var_sq = 1/(bm-am)*sum_temp;
                U_k = sqrt(var_sq).*randn(1,200)+1i*sqrt(var_sq).*randn(1,200);
                u_n_2 = ifft(U_k);
                u_n_2 = ifftshift(abs(u_n_2));
              end
          end
           for k = n*80:(n+1)*80
               s_u_n(1,(n-1)*200+k) = (n+1-k/80)*u_n(k-n*80)+(k/80-n)*u_n_2(k-(n+1)*80);
           end
end
                




