clear
close all
clc

T = 5;
A = 1;
t = linspace(0, T, 50);

s1 = @(t) A*sin(2*pi*t/T).*(t>=0 & t<=T/2);
s2 = @(t) -s1(t-T/2);
 
figure 
subplot(2,1,1)
plot(t, s1(t))
subplot(2,1,2)
plot(t, s2(t))

%%
clear
close all
clc

tic

Eb = 1; % Given value of Eb

SNR = 0:10;
N0 = 10.^(-SNR/10);

Pb_simulated = qfunc(sqrt(10.^(SNR/10))); % Ensure element-wise division by using dot operator

figure
semilogy(SNR, Pb_simulated) % Element-wise division
xlabel('SNR(dB)')
ylabel('Pb (Probability of Bit Error)')
title('Bit Error Rate vs. Eb/N0 (Theoretical Solution for The Same Probabilities)')
grid on

a1 = 1;
a2 = -a1;
gama0 = (a1 + a2)/2;

%bitCount = round(100*1/Pb_simulated(length(SNR))); 
bitCount = 10^8; 
generatedBits = randi([0,1] ,1, bitCount); % generate a1 or a2

ai = generatedBits;
ai(ai==1) = a1;
ai(ai==0) = a2;

Pb_calculated = zeros(1,length(SNR));

n = randn(1,bitCount);

parfor i = 1:length(SNR)
    %bitCount = round(10*1/Pb_simulated(i));
    %generatedBits = randi([0,1], 1, bitCount); % generate a1 or a2
    
    n0 = sqrt(N0(i)).* n; %% calculate noise with variance = N0
    
    z = ai + n0;
    
    s = z > gama0; %% decision
    
    bit_comparison = xor(s, generatedBits); %% calculate how many bits are different
    Pb_calculated(i) = sum(bit_comparison)/bitCount; %% calculate bit error probability by dividing bit count
end


figure
semilogy(SNR, Pb_calculated)
xlabel('SNR(dB)')
ylabel('Pb (Probability of Bit Error)')
title('Bit Error Rate vs. Eb/N0 (Calculated)')
grid on

toc
%%
clear
close all
clc

tic

Eb = 1; % Given value of Eb

SNR = 0:14;
N0 = 10.^(-SNR/10);
a1 = 1;
a2 = -a1;

sigma0 = sqrt(N0*Eb);
gama0 = 0.55*N0;
P1 = 1/4;
P2 = 3/4;

Pb_sim_b = (1 - qfunc((gama0-a1)./sigma0))*P1 + qfunc((gama0-a2)./sigma0)*P2;

figure
semilogy(SNR, Pb_sim_b)
xlabel('SNR(dB)')
ylabel('Pb (Probability of Bit Error)')
title('Bit Error Rate vs. Eb/N0 (Theoretical Solution for Different Probabilities)')
grid on

bitCount = 10^8; 
generatedBits = randsrc(1,bitCount,[1,0; P1,P2]); % Generate bits with different probabilities
ai = generatedBits;
ai(ai==1) = a1;
ai(ai==0) = a2;

Pb_calc_b = zeros(1,length(SNR));

n = randn(1,bitCount);

for i = 1:length(SNR)
    n0 = sqrt(N0(i)).* n; %% calculate noise with variance = N0
    
    z = ai + n0;
    
    s = z > gama0(i); %% decision
    
    bit_comparison = xor(s, generatedBits); %% calculate how many bits are different
    Pb_calc_b(i) = sum(bit_comparison)/bitCount; %% calculate bit error probability by dividing bit count
end

figure
semilogy(SNR, Pb_calc_b)
xlabel('SNR(dB)')
ylabel('Pb (Probability of Bit Error)')
title('Bit Error Rate vs. Eb/N0 (Theoretical Solution for Different Probabilities)')
grid on

toc

%%
