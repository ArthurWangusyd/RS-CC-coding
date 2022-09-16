function code = ccencoder(msg,g1,g2)
msg_l = length(msg);
code = reshape([conv(msg(1:msg_l),g1);conv(msg(1:msg_l),g2)],1,[]);
code = mod(code(1:2*msg_l),2);
end

