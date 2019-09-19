[y,Fs] = audioread('s5.wav');
Tframe = 25;
Tskip = 10;
length_frame = Fs*Tframe/1000;
length_skip = Fs*Tskip/1000;
%frame samples = 200
%shift samples = 80
number_frame = length(y)/length_skip-2;
w = hamming(length_frame);
frames = zeros(number_frame,length_frame);
autocorr = zeros(number_frame,399);
pitch_max = zeros(number_frame);
frames_tf = zeros(number_frame,length_frame);
frames_hamming = zeros(number_frame,length_frame);
sum_phi = 0;
%find pitch estimation
Am = zeros(number_frame,1);
P = zeros(number_frame,1);
errorsum = zeros(number_frame,1);
    for i = 1:number_frame
        frames(i,:) = y(1+(i-1)*length_skip:(i-1)*length_skip+length_frame);
        frames_hamming(i,:) = (frames(i,:).*w').*w';
        %frames_tf(i,:) = fft(frames(i,:).*w');
        autocorr(i,:) = conv(frames_hamming(i,:),fliplr(frames_hamming(i,:)));
        P(i,1) = opt(autocorr(i,:));
        %spectral Envelope

     end

%%
%pitch refinement
fnum = number_frame;
L = length_frame;
Am = zeros(fnum,92);
erro_mat = zeros(fnum,92);
P_new = zeros(size(P));
for n = 1:number_frame
         error = 0;
         original_pn = P(n,1);
         for P_n = (original_pn-2):0.2:(original_pn+2)
             esum = 0;
             Atemp = zeros(1,92);
             for m = 1:P_n-1
                    % preparation
                    Sw = fft(frames(n,:));
                    %Sw = fftshift(Sw);
                    am = ceil((m-1/2)/P_n*L);
                    bm = floor((m+1/2)/P_n*L);
                    % voiced
                    E_v = fft(hamming(L));
                    E_v = fftshift(E_v);
                    E_v = transpose(E_v);
                    temp1 = Sw.*conj(E_v);
                    temp1 = sum(temp1(am:bm));
                    temp2 = abs(E_v).*abs(E_v);
                    temp2 = sum(temp2(am:bm));
                    
                    A_v = temp1/temp2;
                    e_v = (Sw-A_v*E_v).*(Sw-A_v*E_v);
                    e_v = 1/2/pi*sum(e_v(am:bm));
                    % unvoiced
                    E_u = ones(1,200);
                    temp1 = Sw.*E_u;
                    temp1 = sum(temp1(am:bm));
                    temp2 = abs(E_u).*abs(E_u);
                    temp2 = sum(temp2(am:bm));
                    A_u = temp1/temp2;
                    e_u = (Sw-A_u*E_u).*(Sw-A_u*E_u);
                    e_u = 1/2/pi*sum(e_u(am:bm));
                    % compare error
                    if (e_v < e_u)
                        Atemp(1,m) = A_v;
                        em = e_v;
                        erro_mat(n,m) = 2;
                    else
                        Atemp(1,m) = A_u;
                        em = e_u;
                        erro_mat(n,m) = 1;
                    end
                    esum = esum + em;
             end
            if (esum<error)
                error = esum;
                Am(n,:) = Atemp;
                P_new(n,1) = P_n;
            end
         end
end
         %   am = ceil((m-1/2)*(1/P(i,1)*2*pi));
         %   bm = floor((m+1/2)*(1/P(i,1)*2*pi));
            %voiced
         %   Ew = fft(w);
         %   Ew = fftshift(Ew);
         %   Sw = frames(i,:);
         %   nomi = Sw.*conj(Ew);
         %   nomi = sum(nomi(am,bm));
         %   deno = Ew.*Ew;
         %   deno = sum(deno(am,bm));
         %   VAm = nomi/deno;
            %unvocied
            %UEw = ones(200,1);
            %nomi = Sw;
            %nomi = sum(nomi(am,bm));
            %deno = bm-am;
            %UAm = nomi/deno;
            %temp1 = (Sw-UAm).*(Sw-UAm);
            %temp2 = (Sw-VAm*Ew).*(Sw-VAm*Ew);
            %Uerror = 1/2/pi*sum(temp1(am,bm));
            %Verror = 1/2/pi*sum(temp2(am,bm));
            %if(Uerror>Verror)
            %    Am(i,1) = VAm;
            %    error = Verror;
            %else
            %    Am(i,1) = UAm;
            %    error = Uerror;
            %end
            %errorsum(i,1) = errorsum + error;
    
    
    
%%
%spectral envelope







