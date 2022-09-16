%%
clc
clear all
total_nu=20000;
SNR_all = -10:15;
EBN0 = [-5,-4,-3,-2,-1,0,1,2,3,4,5];
message_error_table=zeros(1,length(EBN0));
message_theory_table=zeros(1,length(EBN0));
BER_hard_table=zeros(1,length(EBN0));
BER_soft_table=zeros(1,length(EBN0));

g1=[1 1 1 1 0 0 1];
g2=[1 0 1 1 0 1 1];
L = 7;
tblen = 42;
trellis = cctrellis(L, [171 133]); 

AWGN_Chan=comm.AWGNChannel();
    AWGN_Chan.NoiseMethod='Signal to noise ratio (Eb/No)';
    AWGN_Chan.BitsPerSymbol=log2(2)
    
for snr =EBN0
    AWGN_Chan.EbNo=snr;
   message_error = 0;
   BER_hard = 0;
   BER_soft=0;
   message_theroy_error = 0;
    BER_uncode = 0;
    for nu=1:1:total_nu
       message=randi([0,255],1,239);
        v_x=RS_enconder_255_239(message);
        yu=de2bi(v_x);
        yu=yu';
        for i=1:1:8
            cc_code=ccencoder(yu(i,:),g1,g2);
            sen=cc_code.*2-1;
            %sen = pskmod(cc_code,2); %BPSK modulation
       
            rev= AWGN_Chan(sen);
            rev1 = real(rev)>0;
            %rev1 = pskdemod(rev,2); %demodulation
            
        qwe=quantiz(real(rev),[-1.6,-1.38,-1.15,-0.92,-0.69,-0.46,-0.23,0,0.23,0.46,0.69,0.92,1.15,1.38,1.6]);
            
            r2_x(i,:)= vitdec(qwe,trellis,tblen,'term','soft',4);
        
        r1_x(i,:)= vitdec(rev1,trellis,tblen,'term','hard');

    
       % r2_x(i,:)=vitdec(real(rev_soft),trellis,tblen,'cont','soft');
        end
        r_x_hard=bi2de(r1_x');
        r_x_soft=bi2de(r2_x');
        correct_message_hard=RSdecoder_255_239(r_x_hard');
        correct_message_soft=RSdecoder_255_239(r_x_soft');
        %code_error_number=239-sum(correct_message_hard==message);
        %message_error=message_error+code_error_number;
        BER_hard = BER_hard + sum(sum(abs(de2bi(message)-de2bi(correct_message_hard))));
        BER_soft=BER_soft + sum(sum(abs(de2bi(message)-de2bi(correct_message_soft))));
        %BER_uncode = BER_uncode + sum(sum(abs(de2bi(message)-de2bi(r_uncode))));
    end
   
    
    
    BER_hard_table(snr==EBN0) = BER_hard;
    BER_soft_table(snr==EBN0)=BER_soft;
    fprintf('%d\n',snr);
end
BER_hard_1 = BER_hard_table/(8*239*total_nu);
BER_soft_1=BER_soft_table/(8*239*total_nu);


%%
close all



semilogy(EBN0,BER_hard_1,'^-','LineWidth',2);
hold on
semilogy(EBN0,BER_soft_1,'^-','LineWidth',2);


legend('hard decision Coded','soft decision Coded');
xlabel('EBN0');
ylabel('Bit Error Rate');
title('Bit Error Rate of hard and soft decoder');
