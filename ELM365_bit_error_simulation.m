clear
close all
clc

T = 1;
A = 2; % A and T for Eb = 1
a1 = (A^2*T)/4;
a2 = -a1;
E1 = (A^2*T)/4; % Energy of 1 Bit
E2 = (A^2*T)/4; % Energy of 0 Bit

t = linspace(0, T, 50);

s1 = @(t) A*sin(2*pi*t/T).*(t>=0 & t<=T/2);
s2 = @(t) -s1(t-T/2);
 
figure 
subplot(2,1,1)
plot(t, s1(t))
title("s1(t)")
xlabel("t")
subplot(2,1,2)
plot(t, s2(t))
title("s2(t)")
xlabel("t")

%%
close all
tic

SNR = 0:15;
Pb_simulated = qfunc(sqrt(10.^(SNR/10))); % Ensure element-wise division by using dot operator

figure
semilogy(SNR, Pb_simulated) % Element-wise division
xlabel('SNR(dB)')
ylabel('Pb (Probability of Bit Error)')
title('Bit Error Rate vs. Eb/N0 (Theoretical Solution for The Same Probabilities)')
grid on

P1 = 1/2;
P2 = 1/2;
Eb = E1*P1 + E2*P2;
gama0 = (a1 + a2)/2;
N0 = 10.^(-SNR/10);

bitCount = 10^8; 
generatedBits = randi([0,1] ,1, bitCount); % generate a1 or a2
ai = generatedBits;
ai(ai==1) = a1;
ai(ai==0) = a2;

Pb_calculated = zeros(1,length(SNR));
n = randn(1,bitCount);

parfor i = 1:length(SNR)
    n0 = sqrt(N0(i)).* n; %% calculate noise with variance = N0
    
    z = ai + n0;
    
    s = z > gama0; %% decision 
    
    bit_comparison = xor(s, generatedBits); %% calculate how many bits are different
    Pb_calculated(i) = sum(bit_comparison)/bitCount; %% calculate bit error probability by dividing bit count
end


figure
semilogy(SNR, Pb_calculated, 'b', 'LineWidth', 2)
hold on
semilogy(SNR, Pb_simulated, 'r--', 'LineWidth', 2) 
hold off
xlabel('SNR (dB)')
ylabel('Pb (Probability of Bit Error)')
title('Bit Error Rate vs. Eb/N0')
grid on
legend('Calculated BER', 'Simulated BER')


toc
%%
close all

tic

SNR = 0:15;
N0 = 10.^(-SNR/10);

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
hold on
semilogy(SNR, Pb_sim_b)
hold off
xlabel('SNR (dB)')
ylabel('Pb (Probability of Bit Error)')
title('Bit Error Rate vs. Eb/N0')
grid on
legend('Calculated BER', 'Simulated BER')

toc

%%


figure
semilogy(SNR, Pb_calc_b)
hold on
semilogy(SNR, Pb_calculated)
hold off
xlabel('SNR (dB)')
ylabel('Pb (Probability of Bit Error)')
title('Bit Error Rate vs. Eb/N0')
grid on
legend('Different Probabilities', 'Same Probabilities')



