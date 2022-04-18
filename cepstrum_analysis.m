% Introduction to Speech Processing Project -- Shayan Daneshvar - 9726523
% Calculating Pitch Frequency & Formants Using Cepstrum Transform

clear;
clc;

FL_sec = 0.025; % frame length in seconds
FSL_sec = 0.01; % frame shift length in seconds
lifter_cutoff = 20; % n in liftering of cepstrum

[data,fs]=audioread('a.wav'); % reading the audio

frame_length = FL_sec * fs;
frame_shift_length = FSL_sec * fs;

fn=(length(data)- frame_length)/frame_shift_length + 1; % number of frames
window=hamming(frame_length); % hamming window

pitch_freqs = zeros(1, fn); % pitch frequencies of frames
formants = zeros(4, fn); % 4 formants  of frames
energy_of_frames=zeros(1,fn);% energy of frames


for i=1:fn
    frame = data((i-1)* frame_shift_length + 1:(i-1)* frame_shift_length + frame_length);  %  framing speech
    frame = frame.*window;  %  windowing frame using hamming window
    frame = frame - sum(frame)/length(frame); % removing DC component from frame
    
    quefrency= ifft(log(abs(fft(frame)))); %% cepstrum transform
    
    %%%%% calculating pitch frequncy
    high_liftered = high_time_lifter(quefrency, lifter_cutoff); % high time liftering
    [pks,locs] = findpeaks(high_liftered); % peak picking
    [~,index] = max(pks);
    peak_index = locs(index); %% finding indexes of peaks
    
    if size(peak_index)> 0 % sometimes it is empty
        pitch_freqs(i)= fs/peak_index; %% calculating pitch frequency of frame
    end
     
    %%%%%% calculating formants
    
    low_liftered = low_time_lifter(quefrency, lifter_cutoff); % low-time liftering
    
    low_liftered_dft = abs(log(fft(low_liftered))); % calculating dft to calculate formants
    
    [~,locs] = findpeaks(low_liftered_dft); % peak picking
    
    for formant_no=1:4
        if length(locs)>= formant_no %% returned locs might be less than 4
            formants(formant_no,i) = locs(formant_no) * fs/fn;
        end
    end
    % calculating energy of each frame
    energy_of_frames(i) = sum(frame.*frame)/ length(frame); % calculating frame's energy
end

max_energy = max(energy_of_frames); % finding max energy of a frame
max_energy_frame_index = 0;
for i=1:fn %% finding index of frame with max energy
    if max_energy == energy_of_frames(i)
        max_energy_frame_index = i;
        break
    end    
end  

% recalculating these values for ploting one frame
frame = data(max_energy_frame_index * frame_shift_length + 1:max_energy_frame_index * frame_shift_length + frame_length);
frame = frame.*window;  %  windowing frame using hamming window
frame = frame - sum(frame)/length(frame); % removing DC component from frame
quefrency= ifft(log(abs(fft(frame))));
high_liftered = high_time_lifter(quefrency, lifter_cutoff);
low_liftered = low_time_lifter(quefrency, lifter_cutoff);    
low_liftered_dft = abs(log(fft(low_liftered)));

% ploting cepstrum, high & low liftered , also dft of low time liftered
figure('NumberTitle', 'off', 'Name', 'Frame With Most Energy Information');
subplot(2,2,1);
plot(quefrency); 
title('Cepstrum(Quefrency Domain)');
subplot(2,2,2);
plot(high_liftered); 
title('High-Time-Liftered');
subplot(2,2,3);
plot(low_liftered); 
title('Low-Time-Liftered');
subplot(2,2,4);
plot(low_liftered_dft); 
title('DFT of Low-Liftered');

% Going to Calculate Average pitch and formant between frames with
% acceptable energy, so we get a more reliable result.
avg_pitch_freq = 0;
avg_formant1 = 0;
avg_formant2 = 0;
avg_formant3 = 0;
avg_formant4 = 0;
count_f1 = 0;
count_f2 = 0;
count_f3 = 0;
count_f4 = 0;
count_pitch = 0;

% neglecting frames with energy lower than half the maximum
for i=1:fn
    if energy_of_frames(i) > max_energy/2
        if pitch_freqs(i) > 0
            avg_pitch_freq = avg_pitch_freq + pitch_freqs(i); %% calculating average of pitch frequencies and ignoring 0's due to size(peak_index)> 0
            count_pitch = count_pitch + 1;
        end
        if formants(1,i) > 0
           avg_formant1 = avg_formant1 + formants(1,i); %% calculating average of f1 and ignoring 0's due to length(locs)>= formant_no
           count_f1 = count_f1 + 1;
        end
        if formants(2,i) > 0
           avg_formant2 = avg_formant2 + formants(2,i); %% calculating average of f2 and ignoring 0's due to length(locs)>= formant_no
           count_f2 = count_f2 + 1;
        end
        if formants(3,i) > 0
           avg_formant3 = avg_formant3 + formants(3,i); %% calculating average of f3 and ignoring 0's due to length(locs)>= formant_no
           count_f3 = count_f3 + 1;
        end
        if formants(4,i) > 0
           avg_formant4 = avg_formant4 + formants(4,i); %% calculating average of f4 and ignoring 0's due to length(locs)>= formant_no
           count_f4 = count_f4 + 1;
        end
    end
end

avg_pitch_freq = avg_pitch_freq/count_pitch; %% final result of pitch frequency
avg_formant1 = avg_formant1 / count_f1; %% final result of formants 1,2,3,4
avg_formant2 = avg_formant2 / count_f2;
avg_formant3 = avg_formant3 / count_f3;
avg_formant4 = avg_formant4 / count_f4;

% ploting pitch frequencies
figure('NumberTitle', 'off', 'Name', 'Pitch Frequncies of Frames');
plot(pitch_freqs); 
title('Pitch Frequency');
% ploting formants
figure('NumberTitle', 'off', 'Name', 'Formants of Frames');
hold on
plot(formants(1,:),'r');
plot(formants(2,:),'g');
plot(formants(3,:),'b');
plot(formants(4,:),'c');
title('Formants');

% printing final results (Pitch & Formant)
fprintf("Results: \n Pitch Frequncy %f \n", avg_pitch_freq);
fprintf(" Formants: F1= %f, F2= %f, F3= %f, F4= %f \n", avg_formant1, avg_formant2, avg_formant3, avg_formant4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Functions %%%%%%%%%%%%%

function output = high_time_lifter(frame, n)
    high = zeros(1,length(frame));
    for i=1:length(frame)
        if i >= n
            high(i)=frame(i);
        end    
    end
    output= high;
end

function output = low_time_lifter(frame, n)
    low = zeros(1,length(frame));
    for i=1:length(frame)
        if i < n
            low(i)=frame(i);
        end    
    end
    output= low;
end    
