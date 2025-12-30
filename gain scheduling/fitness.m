function fit=fitness(chr)
global datatrain
global membership
    for i=1:length(datatrain)
                
        c0=membership(i,1)*chr(1)+membership(i,2)*chr(2)+membership(i,3)*chr(3)+membership(i,4)*chr(4)+membership(i,5)*chr(5)+membership(i,6)*chr(6)+membership(i,7)*chr(7)+membership(i,8)*chr(8)+membership(i,9)*chr(9)+membership(i,10)*chr(10)+membership(i,11)*chr(11)+membership(i,12)*chr(12)+membership(i,13)*chr(13)+membership(i,14)*chr(14)+membership(i,15)*chr(15);
        c1=membership(i,1)*chr(16)+membership(i,2)*chr(17)+membership(i,3)*chr(18)+membership(i,4)*chr(19)+membership(i,5)*chr(20)+membership(i,6)*chr(21)+membership(i,7)*chr(22)+membership(i,8)*chr(23)+membership(i,9)*chr(24)+membership(i,10)*chr(25)+membership(i,11)*chr(26)+membership(i,12)*chr(27)+membership(i,13)*chr(28)+membership(i,14)*chr(29)+membership(i,15)*chr(30);
        c2=membership(i,1)*chr(31)+membership(i,2)*chr(32)+membership(i,3)*chr(33)+membership(i,4)*chr(34)+membership(i,5)*chr(35)+membership(i,6)*chr(36)+membership(i,7)*chr(37)+membership(i,8)*chr(38)+membership(i,9)*chr(39)+membership(i,10)*chr(40)+membership(i,11)*chr(41)+membership(i,12)*chr(42)+membership(i,13)*chr(43)+membership(i,14)*chr(44)+membership(i,15)*chr(45);
        G=membership(i,1)*chr(46)+membership(i,2)*chr(47)+membership(i,3)*chr(48)+membership(i,4)*chr(49)+membership(i,5)*chr(50)+membership(i,6)*chr(51)+membership(i,7)*chr(52)+membership(i,8)*chr(53)+membership(i,9)*chr(54)+membership(i,10)*chr(55)+membership(i,11)*chr(56)+membership(i,12)*chr(57)+membership(i,13)*chr(58)+membership(i,14)*chr(59)+membership(i,15)*chr(60);
        
        pred(i,1) = c0+5*log10(datatrain(i,4))+2*log10(datatrain(i,6))+c1*log10(datatrain(i,1)*datatrain(i,6)/datatrain(i,4))+c2*datatrain(i,1)*datatrain(i,6)/datatrain(i,4)+G ;
    end
    errors =abs( (pred - datatrain(:,end)));
    %w=1/(1+(datatrain(:,1).*datatrain(:,6)./datatrain(:,4)).^2);
    %rmse=sum(w.*errors);
    beta = 5;
    
    m = max(errors);
    J2 = m + (1/beta)*log(sum(exp(beta*(errors - m))));
q = 0.95;
J3 = quantile(errors, q);

fit = sum(errors) + 50*J3;


   %fit = sum(errors)+100*max(errors);
   %fit=sum(errors)+J2;
end