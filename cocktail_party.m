clc
close all
%Declaring the audio files
files = {'chirp.mat'
 'gong.mat'};

% Reading/Loading the audio file
S = zeros(10000,2);
for i = 1:2
 test = load(files{i});
 y = test.y(1:10000,1);
 Fs = test.Fs;
 S(:,i) = y;
end

%Mixing the audio files to form a mixed signal
rng default % For reproducibility
mixdata(:,1) = S(:,1)*2 + S(:,2)*1 +5;
mixdata(:,2) = S(:,1)*1 + S(:,2)*2 +5;

%Plotting the mixed signal figure
for i = 1:2
 subplot(2,1,i)
 plot(S(:,i))
 title(['Sound ',num2str(i)])
 subplot(2,1,i)
 plot(mixdata(:,i))
 title(['Mix',num2str(i)])
end


%Extracting the original signal
mixdata = prewhiten(mixdata);
q = 2;
Mdl = rica(mixdata,q,'NonGaussianityIndicator',ones(2,1));
unmixed = transform(Mdl,mixdata);

%Plotting the unmixed signals
figure
for i = 1:2
 subplot(2,1,i)
 plot(S(:,i))
 title(['Sound ',num2str(i)])
 subplot(2,1,i)
 plot(unmixed(:,i))
 title(['Unmixed signal ',num2str(i)])
end

%normalising the unmixed signals
unmixed = unmixed(:,[1:2]);
for i = 1:2
 unmixed(:,i) = unmixed(:,i)/norm(unmixed(:,i))*norm(S(:,i));
end

%Masking
%Assign the seperated signals
mSpeech=unmixed(:,1);
fSpeech=unmixed(:,2);

%Amplitude Adjustment
ampAdj = max(abs([mSpeech;fSpeech;]));
mSpeech = mSpeech/ampAdj;
fSpeech = fSpeech/ampAdj;

% %Assigning the mixed data
mix = mixdata(:,1);

% mix = mSpeech + fSpeech;
mix = mix ./ max(abs(mix));

%Plotting the unmixed and the mixed signals before masking
t = (0:numel(mix)-1)*(1/Fs);
figure()
subplot(3,1,1)
plot(mSpeech)
title("input 1 to binary mask")
grid on
subplot(3,1,2)
plot(fSpeech)
title("input 2 to binary mask")
grid on
subplot(3,1,3)
plot(mix)
title("mixed signal")
xlabel("Time (s)")
grid on

sound(mix);

%Declaring variables for STFT
WindowLength = 128;
FFTLength = 128;
OverlapLength = 96;
win = hann(WindowLength,"periodic");

%Plotting the unmixed and the mixed signal in frequency domain
figure()
subplot(3,1,1)
stft(mSpeech,Fs,'Window',win,'OverlapLength',OverlapLength,'FFTLength',FFTLength)
title("input 1 to binary mask (in frequency domain)")
subplot(3,1,2)
stft(fSpeech,Fs,'Window',win,'OverlapLength',OverlapLength,'FFTLength',FFTLength)
title("input 2 to binary mask (in frequency domain)")
subplot(3,1,3)
stft(mix,Fs,'Window',win,'OverlapLength',OverlapLength,'FFTLength',FFTLength)
title("Mixed audio (in frequency domain)")

%Converting the signals to frequency domain
P_M = stft(mSpeech,'Window',win,'OverlapLength',OverlapLength,'FFTLength',FFTLength);

P_F = stft(fSpeech,'Window',win,'OverlapLength',OverlapLength,'FFTLength',FFTLength);

[P_mix,F] = stft(mix,'Window',win,'OverlapLength',OverlapLength,'FFTLength',FFTLength);

%calculating the binary mask
binaryMask = abs(P_M) >= 2*abs(P_F);

%Applying binary mask to the mixed signal to produce seperate signal
P_M_Hard = P_mix .* binaryMask;
P_F_Hard = P_mix .* (1-binaryMask);

%converting the seperated signal back to time domain
mSpeech_Hard = istft(P_M_Hard ,'Window',win,'OverlapLength',OverlapLength,'FFTLength',FFTLength);
fSpeech_Hard = istft(P_F_Hard,'Window',win,'OverlapLength',OverlapLength,'FFTLength',FFTLength);

%comparing binary mask with original signal
figure()
subplot(2,1,1)
plot(S(:,1))
title("Original Audio 1 Speech")
grid on
mSpeech_Hard = mSpeech_Hard*1.5
subplot(2,1,2)
plot(mSpeech_Hard)
ylim([-1 1])
xlabel("Time (s)")
title("Seperated Audio 1 using binary mask")
grid on
figure()
subplot(2,1,1)
plot(S(:,2))
ylim([-1 1])
title("Original Audio 2 Speech")
grid on
fSpeech_Hard=fSpeech_Hard*2
subplot(2,1,2)
plot(fSpeech_Hard)
ylim([-1 1])
title("Seperated Audio 2 using binary mask")
xlabel("Time (s)")
grid on

%calculating soft mask
softMask = abs(P_M) ./ (abs(P_F) + abs(P_M));

%applying soft mask to the mixed signal to obtain seperated signal
P_M_Soft = P_mix .* softMask;
P_F_Soft = P_mix .* (1-softMask);

%converting the seperated signals back to time domain
mSpeech_Soft = istft(P_M_Soft, 'Window', win, 'OverlapLength', OverlapLength,'FFTLength', FFTLength);
fSpeech_Soft = istft(P_F_Soft, 'Window', win, 'OverlapLength', OverlapLength,'FFTLength', FFTLength);

%comparison of soft mask graphically
figure()
subplot(2,1,1)
plot(S(:,1))
title("Original Audio 1 Speech")
grid on
mSpeech_Soft = mSpeech_Soft*1.5
subplot(2,1,2)
plot(mSpeech_Soft)
ylim([-1 1])
xlabel("Time (s)")
title("Seperated Audio 1 using soft mask")
grid on
figure()
subplot(2,1,1)
plot(S(:,2))
ylim([-1 1])
title("Original Audio 2 Speech")
grid on
fSpeech_Soft = fSpeech_Soft*2
subplot(2,1,2)
plot(fSpeech_Soft)
ylim([-1 1])
title("Seperated Audio 2 using soft mask")
xlabel("Time (s)")
grid on

%Playing Sounds That We Have Seperated
sound(fSpeech_Hard);
sound(mSpeech_Hard);
sound(fSpeech_Soft);
sound(mSpeech_Soft);