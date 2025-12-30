function fit=fitness1st(chr)
global datatrain
global membership
    for i=1:length(datatrain)
        z1 = log10(datatrain(i,1)*datatrain(i,6)/datatrain(i,4));
        z2 = (datatrain(i,1)*datatrain(i,6)/datatrain(i,4));
        
        y_fuzzy = 0;
        
            for k = 1:15
                a0 = chr(3*(k-1)+1);
                a1 = chr(3*(k-1)+2);
                a2 = chr(3*(k-1)+3);
            
                y_fuzzy = y_fuzzy + membership(i,k) * ...
                    (a0 + a1*z1 + a2*z2);
            end
        
        pred(i,1) = ...
            5*log10(datatrain(i,4)) + ...
            2*log10(datatrain(i,6)) + ...
            y_fuzzy;
    end
    errors =abs( (pred - datatrain(:,end)));
    q = 0.95;
    J3 = quantile(errors, q);
    
    fit = sum(errors) + 50*J3;
end