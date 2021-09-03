function [activity_unfiltered, f_max,burst_win,silent_win,prob_fire,amp,period,nu,VMR_act,CV_LFP]=...
    Analysis_fun(N,spike_train,discard,dt,time_length,convolution_bin,LFP,g_syn,sigma,ee,nn)
%% preallocating
[f_max,nu,VMR_act]=deal(zeros(length(g_syn),length(sigma)));

%% fft
%                         spiketrain=spike_train(:,discard/dt:end);
activity=conv(sum(spike_train(:,discard/dt:end)),gausswin(convolution_bin));%/convolution_bin;
activity_unfiltered=activity;

activity2=(activity-(mean(activity))).*(activity>(mean(activity)));
activity2(activity2>0)=1;
% The average line is the baseline for the definition of burst...
% period and silence period. Therefore, we separate the time...
%     intervals above and below this line.
% Burst wins are the time inteval that activity is over average.
% Silent Period are The time interval that activity is below the
% average.
bb=1;
burst_win=0;
silent_win=0;
coh=0;
prob_fire=0;
count=false; %#ok<NASGU>
count_c=false; %#ok<NASGU>
cc=1;
jj=1;
next=0;
[burst_win,silent_win]=deal(-100*ones(1,dt*time_length*0.2));%prealocating 
% With a specific number (-100), preallocating is done to speed up the calculations
% and finally the extra part, by recognizing the specific number, is easily removed. 
% Assuming that the network frequency is ultimately 200 Hz,
% the number of peaks can not exceed this definition.
while jj<length(activity)
    count=1;
    while activity2(jj)==1 && jj<(length(activity)-(convolution_bin)/2)
        % Burst wins are the time inteval that activity is over average.
        % Silent Period are The time interval that activity is below the
        % average.
        
        burst_win(bb)=burst_win(bb)+1;
        coh(bb)=max(coh(bb),activity(jj));
        prob_fire(bb)=prob_fire(bb)+sum(spike_train(:,jj+discard/dt-(convolution_bin)/2-1));
        % firing probability is average number of spike in each burst peiod normalized to the
        % network size.
        
        if bb==1
            true_succession_1=jj;
        end
        jj=jj+1;
        count=0;
        next=1;
    end
    jj=jj+1*count;
    bb=bb+1*next;
    if next==1
        burst_win(bb)=0;
        coh(bb)=0;
        prob_fire(bb)=0;
    end
    next=0;
end
%%
next=0;
jj=1;
while jj<length(activity)
    
    count=1;
    while activity(jj)==0 && jj<length(activity)
        silent_win(cc)=silent_win(cc)+1;
        if cc==1
            true_succession_2=jj;
        end
        jj=jj+1;
        count=0;
        next=1;
    end
    jj=jj+1*count;
    cc=cc+1*next;
    if next
        silent_win(cc)=0;
    end
    next=0;
    
end
%% fft
silent_win=silent_win(silent_win~=-100);% discarding the extra parts remaining from preallocating.
burst_win=burst_win(burst_win~=-100);% discarding the extra parts remaining from preallocating.

if true_succession_1 > true_succession_2% to check the true succession of burst period and silence period
    silent_win=silent_win(2:end);
    
end

T=dt;                    % Sampling frequency
Fs = 1/T;
L=length(activity);
Y = fft(activity);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
conv_power=11;
p=conv(P1,gausswin(conv_power));
[~,f_max_index]=max(p(20:end));

%% results
f_max(ee,nn)=f(f_max_index+20); %Frequency according to maximum poer
[amp,peaks_time]=findpeaks(activity,'MinPeakHeight',mean(activity),'MinPeakDistance',5/dt);
% To find the peaks of population activity as the bursts of network activity.
%Base on the results, the maximum acceptable frequency of
% the network is assumed to be 200 Hz,
% so the distance between the twins should not be less than 5/dt time-steps (5ms).

amp=amp(amp>0.02*N);
% Optionally, it is assumed that a burst is acceptable
% if at least two percent of the network is fired.

period=mean((peaks_time(2:end)-peaks_time(1:end-1))*dt);
CV_LFP=abs(std(LFP)/mean(LFP));
VMR_act(ee,nn)=(mean((activity_unfiltered).^2)-mean(activity_unfiltered)^2)/mean(activity_unfiltered);
nu(ee,nn)=sum(sum(spike_train(:,discard:end)))/(dt*length(spike_train(:,discard:end)));
end
