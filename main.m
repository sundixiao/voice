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
sum_phi = 0;
    for i = 1:number_frame
        frames(i,:) = y(1+(i-1)*length_skip:(i-1)*length_skip+length_frame);
        frames_tf(i,:) = fft(frames(i,:).*w');
        autocorr(i,:) = abs(conv(frames_tf(i,:),fliplr(frames_tf(i,:))));
        for j = 1:200
            for k = 1:399
                if(j*k <= 399)
                    sum_phi = sum_phi + autocorr(i,k*j);
                end
            end
        pitch = sum_phi * j;
            if(pitch > pitch_max(i))
                pitch_max(i) = pitch;
            end
        end
    end

%%
spectrogram(tf,w)
pitch_max = zeros(number_frame);

for j = 1:200
    for k = 1:399
        if(j*k <= 399)
            sum = sum + autocorr(i,k*j);
        end
    end
    pitch = sum * j;
    if(pitch > pitch_max(i))
        pitch_max(i) = pitch;
    end
end
antocorr = xcorr(frames(1,:));
plot(abs(antocorr))

