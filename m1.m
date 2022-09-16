tic
clc
clear 
close

L = 7;
tblen = 42;
trellis = cctrellis(L, [171 133]); %craet  the trellis

%commcnv_plotnextstates(trellis.nextStates); %plot the trellis
%spec = distspec(trellis) %find the dfree and B
SNR = -10:15;%SNR range
len = length(SNR); 
codenum = 10000; % codeword number

g1=[1 1 1 1 0 0 1];
g2=[1 0 1 1 0 1 1];
err = zeros(1,len);
err1 = zeros(1,len); 
err2 = zeros(1,len);%store the error 

codeL = 100; %length of the codeword

for i =1:len
    bit = 0;
    bit1 = 0;
    ratio = 0;
    num = 0;
    num1 = 0;
    bit2 = 0;
    for j = 1:1:codenum
        msg = randi([0 1],1,codeL); %creat the message

        encode=ccencoder(msg,g1,g2); %encode
        

        sen=encode.*2-1; %BPSK modulation
       
        dB = SNR(i) + 10*log10(2)+10*log10(0.5);
        SNR1 = 10^(dB/10);
        noise = randn(1,length(sen));
        noise = noise/sqrt(SNR1); %creat the noise

        
        rev = sen + noise; 
        % AWGN channel
       

        
        

        rev1 = real(rev)>0; %demodulation
        

        decoded = vitdec(rev1,trellis,tblen,'cont','hard'); %hard decision
        num = length(find(msg(1:end-tblen) ~= decoded(tblen+1:end))); %find the different

        bit = bit + num;
       

       


        rev2 = quantiz(real(rev),[-1.6,-1.38,-1.15,-0.92,-0.69,-0.46,-0.23,0,0.23,0.46,0.69,0.92,1.15,1.38,1.6]);
        decode1 = vitdec(rev2,trellis,tblen,'cont','soft',4);
        

        num1 = length(find(msg(1:end-tblen) ~= decode1(tblen+1:end)));
        bit1 = bit1 +num1; % soft decision 

        sen1 = pskmod(msg,2); 

        rev3 = awgn(sen1,SNR(i));
        rev4 = pskdemod(rev3,2);
        bit2 = bit2+ length(find(rev4 ~= msg));
            
    end
    
    err(i) = bit/codeL/codenum;
    err1(i) = bit1/codeL/codenum;
    err2(i) = bit2/codeL/codenum;%calculate the average error
end

figure
semilogy(SNR,err,SNR,err1,SNR,err2,'-o');hold on; %plot the figure 
xlabel('SNR');
ylabel('BER');
legend('hard decision','soft decision','uncoded');
grid on;

toc


