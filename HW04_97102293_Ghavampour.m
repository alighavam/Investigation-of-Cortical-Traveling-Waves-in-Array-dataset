%% Advance Neuro HW04 - Ali Ghavampour - 97102293
clear all; close all; clc;

%% Preprocessing 
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs = 200;
chNum = 48;
trlNum = 823;

% Time Signal
for ch = 1:chNum
    tmp = chan(ch).lfp;
    tmp = mean(tmp,2);
    plot(Time,tmp,'HandleVisibility','off')
    if (ch == 1 || ch == 3)
        plot(Time,tmp,'k','LineWidth',2)
        legend("Bad Channels (1,3)")
    end
    hold on
end
xlim([-1.2 2])
xline(0,'r','LineWidth',2,'HandleVisibility','off');
title("LFP of 48 channels averaged over trials");
xlabel("t(s)")
ylabel("lfp")

figure;
sig = chan(2).lfp;
cleanTrls = Intersect_Clean_Trials;
badTrls = setdiff(1:trlNum,Intersect_Clean_Trials);
for trl = badTrls
    plot(Time,sig(:,trl),'r')
    hold on
end
xlim([-1.2 2])
xline(0,'r','LineWidth',2,'HandleVisibility','off');
for trl = cleanTrls
    plot(Time,sig(:,trl),'k')
    hold on
end
xline(0,'r','LineWidth',2,'HandleVisibility','off');
xlim([-1.2 2])
title("CHANNEL 2 - Black Curves: Clean Trials, Red Curves: Bad Trials");
xlabel("t(s)")
ylabel("lfp")


figure;
x = [2 4 5 6 7 8 9];
for i = 1:length(x)
    ch = x(i);
    sig = chan(ch).lfp;
    plot(Time,mean(sig(:,badTrls),2),'r','LineWidth',1.5)
    hold on
    plot(Time,mean(sig(:,cleanTrls),2),'k','LineWidth',1.5)
    
end
title("Black Curves: avg over clean Trls - Red Curves: avg over bad Trls")
xlabel("t(s)")
ylabel("lfp")
xline(0,'r','LineWidth',2,'HandleVisibility','off');
xlim([-1.2 2])

figure;
for ch = [2,4:chNum]
    tmp = chan(ch).lfp;
    tmp = mean(tmp(:,cleanTrls),2);
    plot(Time,tmp,'HandleVisibility','off')
    hold on
end
xlim([-1.2 2])
xline(0,'r','LineWidth',2,'HandleVisibility','off');
title("LFP of Good channels averaged over clean trials");
xlabel("t(s)")
ylabel("lfp")


figure;
for ch = [1:48]
    sig = chan(ch).lfp;
    sig = mean(sig(:,cleanTrls),2);
    [pxx,f] = pwelch(sig,[],[],[],fs);
    loglog(f,pxx,'k','LineWidth',1.2,'HandleVisibility','off')
    hold on
    if (ch==1 || ch==3)
        loglog(f,pxx,'r','LineWidth',1.2)
        legend("Bad Channels (1,3)")
    end
    hold on
end
title("Power Spectral Density")
xlabel("Frequency (Hz)")
ylabel("Power")


%% Part A & B - Dominant Frequency Using Pwelch
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs = 200;
cleanTrls = Intersect_Clean_Trials;
cleanChs = [2,4:48];

denoisedPow = [];
for ch = cleanChs
    sig = chan(ch).lfp;
    sig = mean(sig(:,cleanTrls),2);
    [pxx,f] = pwelch(sig,[],[0],[],fs);
    x = log10(f);
    y = log10(pxx);
    c = polyfit(x(2:end),y(2:end),1);
    y = c(1)*log10(f) + c(2);
    plot(log10(f),log10(pxx),'k')
    hold on
    plot(log10(f),y,'r')
    hold on
%     plot(log10(f),log10(pxx)-y)
%     hold on
    denoisedPow = [denoisedPow,log10(pxx)-y];
end
fVec = f;
xlim([-0.1 2])
title("Power Spectral Density - Red Lines: linear fit of log power")
xlabel("Log Frequency")
ylabel("Log Power")

figure;
domFreq = [];
for ch = 1:size(denoisedPow,2)
    plot(fVec(:),denoisedPow(:,ch))
    hold on
    max2 = max(denoisedPow(denoisedPow(1:60,ch)<max(denoisedPow(1:60,ch)),ch));
    max1 = max(denoisedPow(1:55,ch));
    ind = find(denoisedPow(:,ch) == max1);
    domFreq = [domFreq,fVec(ind)];
end
xlabel("Frequency")
ylabel("Log Power")
title("Denoised Power Spectral Density of Good Channels")

gp = zeros(size(ChannelPosition));
figure;
for i = 1:size(gp,1)
    for j = 1:size(gp,2)
        ch = ChannelPosition(i,j);
        if (~isnan(ch) && ch~=1 && ch~=3)
            ind = find(cleanChs == ch);
            gp(i,j) = domFreq(ind);
        end
    end
end
% colormap(jet(4))
imagesc(gp)
set(gca,'YDir','normal')
colorbar

figure;
scatter(1:46,domFreq)

%% Dominant Frequency Using FFT
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs = 200;
cleanTrls = Intersect_Clean_Trials;
cleanChs = [2,4:48];

denoisedFFT = [];
for ch = cleanChs
    sig = chan(ch).lfp;
    sig = mean(sig(:,cleanTrls),2);
    [f, mag] = myfft(sig', fs);
    x = log10(f);
    y = log10(mag);
    c = polyfit(x(2:end),y(2:end),1);
    y = c(1)*log10(f) + c(2);
    plot(log10(f),log10(mag),'k')
    hold on
    plot(log10(f),y,'r')
    hold on
    denoisedFFT = [denoisedFFT;log10(mag)-y]; % DomFreqFFT01
%     denoisedFFT = [denoisedFFT; mag-10.^y]; % DomFreqFFT02
end
fVec = f;
xlim([-0.1 2])
title("FFT - Red Lines: linear fit of log fft")
xlabel("Log Frequency")
ylabel("Log fft")

figure;
domFreq = [];
for ch = 1:size(denoisedFFT,1)
    plot(fVec,denoisedFFT(ch,:))
    hold on
    max2 = max(denoisedFFT(ch,denoisedFFT(ch,1:145)<max(denoisedFFT(ch,1:145))));
    max1 = max(denoisedFFT(ch,1:145)); % DomFreqFFT01
%     max1 = max(denoisedFFT(ch,17:145)); % DomFreqFFT02
    ind = find(denoisedFFT(ch,:) == max1);
    domFreq = [domFreq,fVec(ind)];
end
xlabel("Frequency")
ylabel("Log fft")
title("Denoised Power Spectral Density of Good Channels")

gp = zeros(size(ChannelPosition));
figure;
for i = 1:size(gp,1)
    for j = 1:size(gp,2)
        ch = ChannelPosition(i,j);
        if (~isnan(ch) && ch~=1 && ch~=3)
            ind = find(cleanChs == ch);
            gp(i,j) = domFreq(ind);
        end
    end
end
imagesc(gp)
set(gca,'YDir','normal')
colorbar

figure;
scatter(1:46,domFreq)



%% Part C - Spectrogram - Welch Method
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs = 200;
cleanTrls = Intersect_Clean_Trials;
cleanChs = [2,4:48];

tTiles = 21;
inds = round(linspace(1,length(Time),tTiles));
pMat = zeros(10,length(inds)-1);
for ch = cleanChs
    sig = chan(ch).lfp;
    sig = mean(sig(:,cleanTrls),2);
    tVec = [];
    pMatTmp = [];
    for i = 1:length(inds)-1
        tmp = sig(inds(i):inds(i+1));
        [pxx,f] = pwelch(tmp,[],[],[33],fs);
        x = log10(f);
        y = log10(pxx);
        c = polyfit(x(2:end),y(2:end),1);
        y = c(1)*log10(f) + c(2);
        tmpInd = find(f>=50);
        tmpInd = tmpInd(1);
        f = f(1:tmpInd);
        pxx = pxx(1:tmpInd); 
        pMatTmp = [pMatTmp,pxx];
        tVec = [tVec, Time(inds(i))];
    end
    pMat = pMat + pMatTmp;
end
imagesc(tVec*1000,f,pMat/(length(inds)-1));
title("power spectrogram - Welch's method")
xlabel("t(ms)")
ylabel("frequency(Hz)")
colorbar
set(gca,'YDir','normal')
xline(0,'k','LineWidth',2);

%% Spectrogram - Multitaper Method
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs = 200;
cleanTrls = Intersect_Clean_Trials;
cleanChs = [2,4:48];

tTiles = 31;
inds = round(linspace(1,length(Time),tTiles));
pMat = zeros(12,length(inds)-1);
for ch = cleanChs
    sig = chan(ch).lfp;
    sig = mean(sig(:,cleanTrls),2);
    tVec = [];
    pMatTmp = [];
    for i = 1:length(inds)-1
        tmp = sig(inds(i):inds(i+1));
        [pxx,f] = pmtm(tmp,[3],[length(tmp)],fs);
        pMatTmp = [pMatTmp,pxx];
        tVec = [tVec, Time(inds(i))];
    end
    pMat = pMat + pMatTmp;
end
colormap(jet)
imagesc(tVec*1000,f,pMat/(length(inds)-1));
title("power spectrogram - Multitaper method")
xlabel("t(ms)")
ylabel("frequency(Hz)")
colorbar
set(gca,'YDir','normal')
xline(0,'k','LineWidth',2);

%% Channels Spectrogram
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs = 200;
cleanTrls = Intersect_Clean_Trials;
cleanChs = [2,4:48];

denoisedPow = [];
noisedPow = [];
for ch = cleanChs
    sig = chan(ch).lfp;
    sig = mean(sig(:,cleanTrls),2);
    [pxx,f] = pwelch(sig,[],[0],[],fs);
    x = log10(f);
    y = log10(pxx);
    c = polyfit(x(2:end),y(2:end),1);
    y = c(1)*log10(f) + c(2);
    noisedPow = [noisedPow,log10(pxx)];
    denoisedPow = [denoisedPow,log10(pxx)-y];
end

indf = find(f>=45);
indf = indf(1);
f = f(1:indf);
noisedPow = noisedPow(1:indf,:);
denoisedPow = denoisedPow(1:indf,:);

figure
subplot(1,2,1)
colormap(jet)
imagesc(1:46,f,noisedPow)
xlabel("Channels")
ylabel("Frequency")
title("Noisy Power")
colorbar
set(gca,'YDir','normal')
subplot(1,2,2)
colormap(jet)
imagesc(1:46,f,denoisedPow)
xlabel("Channels")
ylabel("Frequency")
title("Denoised Power")
colorbar
set(gca,'YDir','normal')

%% Phase Propagation ===========================================================
clear all; close all; clc;
load('ArrayData.mat')
load('CleanTrials.mat')
fs = 200;
cleanTrls = Intersect_Clean_Trials;
cleanChs = [2,4:48];
filtOrder = 2;

domFreq = load('domFreqFFT01'); % Choose DomFreq Here !!!!!!!!!!!!!!!!!
domFreq = domFreq.domFreq;

filtSig = zeros(length(cleanChs),length(cleanTrls),size(chan(2).lfp,1));
phi = zeros(length(cleanChs),length(cleanTrls),size(chan(2).lfp,1));
tic
for i = 1:length(cleanChs)
    ch = cleanChs(i);
    fc = domFreq(i);
    for j = 1:length(cleanTrls)
        trl = cleanTrls(j);
        tmpSig = chan(ch).lfp;
        tmpSig = tmpSig(:,trl);
        fcut1 = fc - 1;
        fcut2 = fc + 1;
        [b,a] = butter(filtOrder,[fcut1,fcut2]/(fs/2),'bandpass');
        filtTmpSig = filter(b,a,tmpSig);
        filtSig(i,j,:) = filtTmpSig;
        tmpPhi = instPhase(filtTmpSig);
        phi(i,j,:) = tmpPhi;
    end
end
toc

%% Wave Animation Demo - NEED TO RUN "PHASE PROPAGATION" SECTION!
clc
close all
targTrial = randi(length(cleanTrls))
targTrial = 38;
map = zeros(size(ChannelPosition,1),size(ChannelPosition,2),size(chan(2).lfp,1));
for i = 1:size(map,1)
    for j = 1:size(map,2)
        ch = ChannelPosition(i,j);
        if (~isnan(ch) && ch ~= 1 && ch ~= 3)
            ind = find(cleanChs == ch);
            tmpPhi = reshape(phi(ind,targTrial,:),size(chan(2).lfp,1),1);
            tmpCos = cos(tmpPhi);
            map(i,j,:) = tmpCos;
        end
    end
end
map = map(:,2:end,:);

t1 = 201;
t2 = 321;
h = figure('Position',[500 200 1000 500]);
colormap(gray)
pgd = PGD(map);
spd = speed(map,fs);
v = VideoWriter('myVideo.avi');
open(v)
for t = t1:t2
    if ~ishghandle(h)
        break
    end
    subplot(2,5,[1:3,6:8])
    imagesc(map(:,:,t));
    set(gca,'YDir','normal')
    hold on
    [FX,FY] = gradient(map(:,:,t));
    x = 1:9;
    y = 1:5;
    quiver(x,y,FX,FY);
    title(sprintf("t = %dms , Trial %d",round(Time(t)*1000),cleanTrls(targTrial)))
    
    subplot(2,5,9)
    bar(pgd(t))
    title("PGD")
    ylim([0 1])
    
    subplot(2,5,10)
    bar(spd(t))
    title("Speed")
    ylabel("cm/s")
    ylim([0 100])
    
    subplot(2,5,[4,5])
    theta = atan(reshape(FY,1,45) ./ reshape(FX,1,45));
    polarhistogram(theta,11);
    title("Direction Histogram")
    
    frame = getframe(gcf);
    for i = 1:5
        writeVideo(v,frame);
    end
%     pause(0.3)
end
close(v)

figure;
t = round(linspace(241,251,6));
for i = 1:6
    subplot(2,3,i)
    imagesc(map(:,:,t(i)));
    set(gca,'YDir','normal')
    title(sprintf("t = %d , trial %d",round(Time(t(i))*1000),targTrial))
end

%% Direction Histogram - NEED TO RUN "PHASE PROPAGATION" SECTION!
clc
close all
fxVec = [];
fyVec = [];
t1 = 201;
t2 = 321;
tic
for trl = 1:length(cleanTrls(1:20))
    targTrial = cleanTrls(trl);
    map = zeros(size(ChannelPosition,1),size(ChannelPosition,2),t2-t1+1);
    for i = 1:size(map,1)
        for j = 1:size(map,2)
            ch = ChannelPosition(i,j);
            if (~isnan(ch) && ch ~= 1 && ch ~= 3)
                ind = find(cleanChs == ch);
                tmpPhi = reshape(phi(ind,targTrial,:),size(chan(2).lfp,1),1);
                tmpPhi = tmpPhi(t1:t2);
                tmpCos = cos(tmpPhi);
                map(i,j,:) = tmpCos;
            end
        end
    end
    map = map(:,2:end,:);
    [FX,FY] = gradient(map);
    pgd = PGD(map);
    for t = 1:t2-t1+1
%         if (pgd(t)>=0.5)
%             fxVec = [fxVec, reshape(FX(:,:,t),1,45)];
%             fyVec = [fyVec, reshape(FY(:,:,t),1,45)];
%         end
        fxVec = [fxVec, reshape(FX(:,:,t),1,45)];
        fyVec = [fyVec, reshape(FY(:,:,t),1,45)];
    end
end
toc

% Direction Histogram
dirVec = [];
for i = 1:length(fxVec)
    x = fxVec(i);
    y = fyVec(i);
%     if (x>=0) % -180 to 180
%         dirVec(i) = atand(y/x);
%     elseif (x<0 && y>=0)
%         dirVec(i) = atand(y/x) + 180;
%     elseif (x<0 && y<0)
%         dirVec(i) = atand(y/x) - 180;
%     end
%     if (dirVec(i) <= -165)
%         dirVec(i) = -dirVec(i);
%     end
    dirVec(i) = atand(y/x); % -90 to 90
end

figure;
% x = -165:180; % -180 to 180
x = -90:90;
h = hist(dirVec,x);
bar(x,h,'k')
title("Gradient Direction Histogram")
xlabel("Angle")
ylabel("Count")


% Speed Histogram
spdVec = [];
t1 = 201;
t2 = 321;
tic
for trl = 1:length(cleanTrls(1:20))
    targTrial = cleanTrls(trl);
    map = zeros(size(ChannelPosition,1),size(ChannelPosition,2),t2-t1+1);
    for i = 1:size(map,1)
        for j = 1:size(map,2)
            ch = ChannelPosition(i,j);
            if (~isnan(ch) && ch ~= 1 && ch ~= 3)
                ind = find(cleanChs == ch);
                tmpPhi = reshape(phi(ind,targTrial,:),size(chan(2).lfp,1),1);
                tmpPhi = tmpPhi(t1:t2);
                tmpCos = cos(tmpPhi);
                map(i,j,:) = tmpCos;
            end
        end
    end
    map = map(:,2:end,:);
    spd = speed(map,fs);
    spdVec = [spdVec, spd];
end
toc

% Histogram
figure;
x = 0:200;
h = hist(spdVec,x);
bar(x,h,'k')
title("Speed Histogram")
xlabel("speed(cm/s)")
ylabel("Count")



%% Functions
function spd = speed(map,fs)
    spd = zeros(1,size(map,3));
    for t = 1:length(spd)-1
        [fx,fy] = gradient(map(:,:,t));
        [fx2,fy2] = gradient(map(:,:,t+1));
        
        tmpX = mean(fx2 - fx,'all')*fs;
        tmpY = mean(fy2 - fy,'all')*fs;
        num = sqrt(tmpX^2 + tmpY^2);
        
        normTmp = [];
        for i = 1:5
            for j = 1:9
                normTmp = [normTmp, norm([fx(i,j), fy(i,j)])];
            end
        end
        denom = mean(normTmp);
        spd(t) = num/denom;
    end 
end

function pgd = PGD(map)
    pgd = zeros(1,size(map,3));
    for t = 1:length(pgd)
        [fx,fy] = gradient(map(:,:,t));
        avgX = mean(fx,'all');
        avgY = mean(fy,'all');
        num = norm([avgX, avgY]);
        normTmp = [];
        for i = 1:5
            for j = 1:9
                normTmp = [normTmp, norm([fx(i,j), fy(i,j)])];
            end
        end
        denom = mean(normTmp);
        pgd(t) = num/denom;
    end
end

function phase = instPhase(sig)
    x = hilbert(sig);
    Sa = sig + 1j*x;
    phase = angle(Sa);
end

function [f, mag] = myfft(x, fs)
    y = fft(x);
    mag = abs(y);
    mag = fftshift(mag)/fs;
    L = size(x,2);
    f = linspace(-fs/2,fs/2,L);
    ind = find(f >= 0);
    ind = ind(1);
    f = f(ind:end);
    mag = mag(ind:end);
end

