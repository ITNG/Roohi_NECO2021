function[spike_train,synch_index,LFP]=model_simulator(N,time_length,dt,model,tau_decay,delay,p_conn_II,sigma,g_syn,nn,ee,discard)

I_ext=2;% external drive (micro Ampere)
time=0:dt:time_length; %milisecond
warning('off')

%% Constant
%-------- Networks parameters---------
N_E=0; N_I=N-N_E;
not_test=ones(N,1);
not_test(10)=1;

% ---------Neuron temporal parameters----------
tau_m=10;%Membrane time constant
tau_rise_ampa=0.5; %Tau_decay=5; %milisecond
tau_rise_gaba_a=0.5;%milisecond

%% function

if model==1
    % ----------------------1- LIF-------------------
    
    v_rest=-70; v_th=-50;%mv
    R=10;%Kiloo_ohm
    vdot=@(v,I_external,I_synaptic,I_noise) tau_m\...
        (-(v-v_rest)+R*(I_external+I_synaptic+I_noise));
    
elseif model==2
    
    % ----------------------3-EIF--------------------
    tau_m= 10;%ms
    Delta_T=1; v_rest=-70; v_rh=-50; v_th=20;%milivolt
    R=10;% Kiloo_Ohm
    vdot=@(v,I_external,I_synaptic,I_noise) tau_m\(-(v-v_rest)+...
        Delta_T*exp((v-v_rh)/Delta_T)+R*(I_external+I_synaptic+I_noise));
    
end

%% Connectivity
%connection probabilty
EE=1; EI=1; IE=1; II=1*p_conn_II;
var_conn=0; % standard deviation of connectivity inhomogeneity
A=create_A(EE,EI,IE,II,N_E,N_I);   w_connection=reshape(normrnd(1,var_conn,N*N,1),N,N);
A=A.*w_connection;

%% Synapse: base on the exponentialy rise-decay behavior

S_ij=@(t,t_f,tau_rise,tau_decay)...
    (((1-exp(-tau_rise\((t-delay)-t_f))).*exp(-tau_decay\((t-delay)-t_f)).*(sign((t-delay)-t_f)+1)/2));
max_S_ij_Excitatory=max(S_ij(delay:dt:10+delay,0,tau_rise_ampa,tau_decay));
max_S_ij_inhibitory=max(S_ij(delay:dt:10+delay,0,tau_rise_gaba_a,tau_decay));

%% Parameters Pre-definition

I_syn_0=zeros(N,1);%initial value
v_0=normrnd(v_rest+I_ext*R,0.8,N,1);%initial random value of neurons voltage
fired=false(N,1); %#ok<PREALL>
last_fire=zeros(N,1); %#ok<PREALL>
spike_train=false(N,length(time)); %#ok<PREALL>

I_noise=normrnd(0,sigma(nn),N,length(time));
epsilon_syn=g_syn(ee)*ones(N,1)/N;
I_syn=I_syn_0;
v=v_0;
spike_train=false(N,length(time)); fired=false(N,1);  last_fire=zeros(N,1); %#ok<PREALL>
%% The Numeric Calculation

%% preallocating for eachloop
LFP=zeros(1,length(time));
count_to_average=0;
A_N=zeros(1,length(time)-discard/dt);
V_1=zeros(N,1);
V_2=zeros(N,1);

for tt=1:length(time)-1
   
    %% --------Eulerian Method--------
    
    v=dt*vdot(v,I_ext,I_syn,not_test.*I_noise(:,tt))+v;
    fired=v>v_th;
    v(fired)=v_rest;
    spike_train(:,tt)=fired;
    last_fire=max(tt*fired,last_fire);
    t=tt;
    
    % synchrony index parameters
    if tt>(discard/dt)
        A_N(tt-discard/dt)=mean(v);
        V_2=V_2+v.^2;
        V_1=V_1+v;
        count_to_average=count_to_average+1;
        %Delta_N=var(A_N);
    end
    
    
    %% --------synaptic current calculation base on the synaptic profile--------------
    I_syn_inhibitory=20*(-max_S_ij_inhibitory\A(:,N_E+1:N)*((epsilon_syn(N_E+1:N)).*...
        (S_ij(t*dt,last_fire(N_E+1:N)*dt,tau_rise_gaba_a,tau_decay))));% inhibitory synaptic current
    
    I_syn_excitatory=20*(max_S_ij_Excitatory\A(:,1:N_E)*((epsilon_syn(1:N_E).*(~(~last_fire(1:N_E)))).*...
        S_ij(t*dt,last_fire(1:N_E)*dt,tau_rise_ampa,tau_decay)));% excitatory synaptic current.
    % last part is to check whether neuron start firing or not.
    I_syn=I_syn_excitatory+I_syn_inhibitory;%total synaptic current
    LFP(tt)=sum(sum(I_syn)*R);% Average synaptic current assumed as the Local Field Potential (LFP)
end

%% synchrony index calculation
A_N= A_N(1:end-1);
Delta_N=mean(A_N.^2)-(mean(A_N))^2;
Delta=mean((V_2)/count_to_average-((V_1)/count_to_average).^2);
synch_index=Delta_N/Delta;
end
