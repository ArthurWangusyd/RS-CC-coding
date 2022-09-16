
function codeword = RS_enconder_255_239(msg)
    m = 8;
    t=8;
    
    % if use the X^(2t)*c(x), it will left shift 2t
     
    x_2t = [ msg zeros(1,2*t)];
    %  g(x) and convert x_2t to gf
    gen_poly = rsgenpoly(255,239);
    x_2tc = gf(x_2t,m);
    %  b(x) 
    [a,b] = deconv(x_2tc,gen_poly);
    % v(x) = X^(2t)*c(x) + b(x)
    v = x_2tc(1:255) + b(1:255); 
    codeword = double(v.x);
end