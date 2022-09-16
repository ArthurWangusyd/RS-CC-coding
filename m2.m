% this file will plot the relationship between pe and BER
clc
clear 
close
tic

pe = 0:0.001:0.01;
pe1 = 0.02:0.01:0.1;
pe = [pe pe1]; % range of pe
len = length(pe);



L = 7;
tblen = 42;
trellis = cctrellis(L, [171 133]); %craet  the trellis

%commcnv_plotnextstates(trellis.nextStates); %plot the trellis
%spec = distspec(trellis) %find the dfree and B

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

        encode=ccencoder(msg,g1,g2);
        

        
       
        rev = bsc(encode,pe(i));
        msg1 = bsc(msg,pe(i));
        
        
        
        

        decoded = vitdec(rev,trellis,tblen,'cont','hard'); %hard decision
        num = length(find(msg(1:end-tblen) ~= decoded(tblen+1:end))); %find the different

        bit = bit + num;

        rev2 = quantiz(rev,[0,1]);
        decode1 = vitdec(rev2,trellis,tblen,'cont','soft',1);
        
        
        num1 = length(find(msg(1:end-tblen) ~= decode1(tblen+1:end)));
        bit1 = bit1 +num1; % soft decision 

        bit2 = bit2+ length(find(msg1 ~= msg));
            
    end
    
    err(i) = bit/codeL/codenum;
    err1(i) = bit1/codeL/codenum;
    err2(i) = bit2/codeL/codenum;%calculate the average error
end
figure
%%
semilogy(pe,err)
hold on
semilogy(pe,err1,"--")
hold on
semilogy(pe,err2)
hold on
xlabel('Pe');
ylabel('BER');
legend('hard decision','soft decision','uncoded');
grid on;



toc