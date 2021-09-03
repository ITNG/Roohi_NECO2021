clear
close all
rng('shuffle')
maxNumCompThreads(3)% set your appropriate number of threads based on your computer or cluster 
% The codes, contain the "Main" part where you can specify
% the values of the simulation parameters.
% In function "Create_A", the connection matrix 
%is made based on the specified values.
%In the "Model_simulator", numerical calculations are 
%performed and the output is the spike-train matrix.
%In the "Analaysis_fun", the analyzes are performed on the 
%results and all the required simulation outputs are extracted.
% To repeat the simulations, you can set the parameter in the 
%"Main" and run it in \MATLAB (2019a or later versions).

%% ------------------Parameters------------
Tau_decay=1:2:5;% synaptic decay time (ms)
DELAY=[1.5,2,2.5,3,3.5,4,4.5,5];% transmision delay (ms)
time_length=20000;% simulation time length
dt=0.1;  %Time step of simulation (ms)
model=1;% Select the model:1-LIF, 2-EIF
N=1024;% Network size
conv_bin=5/dt; % time length of Gaussian window of convolution
discard=300;% time interval to discard the initial condition
p_conn_II=50/100;% connection probabilty
sigma=linspace(0.5,25,10);% standard deviation of the  Gausian noise
g_syn=p_conn_II\linspace(1.5,12.5,10);% synaptic coupling strength

sigma=sigma(4);
g_syn=g_syn(5);
% If you want to do the simulation for the whole parameter space, comment
% the two lines above.
%% -----------------loops for decay time and delay------------------------
for td=2%1:length(Tau_decay)
    tau_decay=Tau_decay(td);
    
    for dl=3%1:length(DELAY)
        delay=DELAY(dl);
        
        %% ---------------loops for noise strength and coupling strnegth
        for ee=1:length(g_syn)
            for nn=1:length(sigma)
                
                %% --------------------Simulation----------------
                [spike_train,synch_index,LFP]=...
                    model_simulator(N,time_length,dt,model,tau_decay,delay,p_conn_II,sigma,g_syn,nn,ee,discard);
                
                %% -------------------Analysis-------------------
                [activity_unfiltered, f_max,burst_win,silent_win,prob_fire,amp,period,nu,VMR_act,CV_LFP]=...
                    Analysis_fun(N,spike_train,discard,dt,time_length,conv_bin,LFP,g_syn,sigma,ee,nn);
                
                %% ------------------------Save-------------------
                
                save(['results_delay_test',num2str(dl),'_tau',num2str(td)],'f_max',...
                    'period','nu','LFP','activity_unfiltered','CV_LFP','synch_index','VMR_act',...
                    'sigma','g_syn','burst_win','silent_win','prob_fire','amp')
                
                
            end
            
        end
    end
end

