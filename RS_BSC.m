%%
clc
clear all
total_nu=500000;
p_all = logspace(-4,-1,8);
message_error_table=zeros(1,length(p_all));
message_theory_table=zeros(1,length(p_all));
BER_coded_table=zeros(1,length(p_all));
BER_uncoded_table=zeros(1,length(p_all));

for p =p_all(1:8)
   message_error = 0;
   BER = 0;
   message_theroy_error = 0;
    BER_uncode = 0;
    for nu=1:1:total_nu
       message=randi([0,255],1,239);
        v_x=RS_enconder_255_239(message);
        yu=bi2de(bsc(de2bi(v_x),p));
        r_x=(yu)';
        correct_message=RSdecoder_255_239(r_x);
        r_uncode=bi2de(bsc(de2bi(message),p))';
        
        code_error_number=239-sum(correct_message==message);
        message_error=message_error+code_error_number;
        BER = BER + sum(sum(abs(de2bi(message)-de2bi(correct_message))));
        BER_uncode = BER_uncode + sum(sum(abs(de2bi(message)-de2bi(r_uncode))));
    end
    p_byte = 1-(1-p)^8;
    for j = 0:8
        message_theroy_error = message_theroy_error+ nchoosek(255-j,j)*(p_byte)^j*(1-p_byte)^(255-j);
    end
    message_theory_table(p==p_all) = 1 - message_theroy_error;
    message_error_table(p==p_all) = message_error;
    BER_coded_table(p==p_all) = BER;
    BER_uncoded_table(p==p_all) = BER_uncode;
    fprintf('%d\n',p);
end

MEP = message_error_table/(239*total_nu);
BER = BER_coded_table/(8*239*total_nu);
BER_uncode = BER_uncoded_table/(8*239*total_nu);

%%
close all
figure (1)
semilogy(p_all,MEP,'^-','LineWidth',2);
hold on
semilogy(p_all,message_theory_table,'s-','LineWidth',2);
legend('Simulation','Theory');
xlabel('BSC Probability p');
ylabel('Message Error Probability');
title('Message Error Probability of simulation vs actual');

figure (2)
semilogy(p_all,BER,'^-','LineWidth',2);
hold on
semilogy(p_all,BER_uncode,'s-','LineWidth',2);
legend('Coded','Uncoded');
xlabel('BSC Probability p');
ylabel('Bit Error Rate');
title('Bit Error Rate of coded and uncoded');
