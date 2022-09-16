

function message = RSdecoder_255_239(codeword)
    
 
    
    % extend r to 255
    re=codeword;
    codeword = gf(re,8);
    % s(x) be calculated, since t=8, s_2t should have 16 num
    alpha = gf(2,8);
    s_x = polyval(codeword,alpha.^(1:16));
    
    % B-MAlgorithm
    sigma_miu = [gf(zeros(1,8),8) 1];
    %initila result==1
    d_miu = gf(1,8);
    sigma_plus = sigma_miu; % sigma^0_(x)
   
    miu_l = [-1;-1];
    sigma_list = sigma_miu;
    d_list = d_miu;
    
    % check s(x) is whether all zeros
    if s_x == gf(zeros(1,255-239),8)
        message = double(codeword.x(1:255-255+239));
    else
        for miu = 0:255-239-1
            sigma_miu = sigma_plus;
            %  I_miu is the degree of sigma^miu_(x)
            yuyuyu=min(find(sigma_miu ~= 0));
            I_miu = 8-yuyuyu+1;
                 
  
            
            % first to calculate d_miu
            d_miu = s_x(miu+1);
            if I_miu >= 1
                for j = 1:I_miu
                    d_miu = d_miu + s_x(miu+1-j)*sigma_miu(end-j);
                end
            end
            if d_miu ~= 0
                % find rho
                rho_idx = find(miu_l(1,:)==max(miu_l(1,:)));
                rho = miu_l(2,rho_idx(end));
                
                miu_l = [miu_l [miu-I_miu;miu]];
                
                sigma_rho = sigma_list(rho+2,:);
                d_rho = d_list(rho+2);
                % caculate the sigma_(miu+1)
                add = d_miu/d_rho*[sigma_rho zeros(1,miu-rho)];
                sigma_plus = sigma_miu + add(end-8:end);
                % the max error that can correct is (N-K)/2
                if I_miu > 8
                    message = double(codeword.x(255-255+1:255-255+239));
                    return
                end
            end
            
            sigma_list = [sigma_list;sigma_miu];
            d_list = [d_list d_miu];
        end
        sigma = sigma_plus; % final sigma(x)
        % calculate the roots of sigma(x),to find error location which is 1/beta
        beta_1 = (roots(sigma))';
        if size(beta_1,1) == 0
            message = double(codeword.x(255-255+1:255-255+239));
            return
        end
        GF = gftuple((-1:2^8-2)',8);
        check_table = zeros(1,2^8);
        for k = 1:2^8
            check_table(k) = bi2de(GF(k,:));
        end
        % obtain the error location
        error_index = zeros(1,size(beta_1,2));
        temp = double(beta_1.x);
        for a = 1:size(beta_1,2)
            if temp(a) == 1
                error_index(a) = 255;
            else
                error_index(a) = find(check_table==temp(a))-2;
            end
        end
        % remember, s(x) from left to right and sigma from right to left
        Z = [s_x(8)+sigma(8)*s_x(7)+sigma(7)*s_x(6)+sigma(6)*s_x(5)+sigma(5)*s_x(4)+sigma(4)*s_x(3)+sigma(3)*s_x(2)+sigma(2)*s_x(1)+sigma(1)
             s_x(7)+sigma(8)*s_x(6)+sigma(7)*s_x(5)+sigma(6)*s_x(4)+sigma(5)*s_x(3)+sigma(4)*s_x(2)+sigma(3)*s_x(1)+sigma(2)
             s_x(6)+sigma(8)*s_x(5)+sigma(7)*s_x(4)+sigma(6)*s_x(3)+sigma(5)*s_x(2)+sigma(4)*s_x(1)+sigma(3)
             s_x(5)+sigma(8)*s_x(4)+sigma(7)*s_x(3)+sigma(6)*s_x(2)+sigma(5)*s_x(1)+sigma(4)
             s_x(4)+sigma(8)*s_x(3)+sigma(7)*s_x(2)+sigma(6)*s_x(1)+sigma(5)
             s_x(3)+sigma(8)*s_x(2)+sigma(7)*s_x(1)+sigma(6)
             s_x(2)+sigma(8)*s_x(1)+sigma(7)
             s_x(1)+sigma(8) 
             1];
        %find the value of e
        error_value = gf(zeros(1,size(beta_1,2)),8);
        for x = 1:size(beta_1,2)
            error_value(x) = polyval(Z,beta_1(x));
            for y = 1:size(beta_1,2)
                if beta_1(x) ~= beta_1(y)
                    error_value(x) = error_value(x) / (beta_1(x)/beta_1(y) + 1);
                end
            end
        end
        % trans e to  e(x)
        err = gf(zeros(size(codeword)),8);
        err(error_index) = error_value;
        % find the v(x) 
        v_x = codeword + err;
        % check the s_x
        if polyval(v_x,alpha.^(1:255-239)) == gf(zeros(1,255-239),8)
            message = double(v_x.x((255-255+1:255-255+239)));
        else
            message = double(codeword.x(255-255+1:255-255+239));
        end
    end
    
end