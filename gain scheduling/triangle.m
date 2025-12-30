function mu = triangle(x,le,ce,re)    %TRIANGLE MEMBERSHIP FUNCTION
    if x >= le & x < ce
       mu = (x - le)/(ce - le);
    elseif x >= ce & x <= re
       mu = (re - x)/(re - ce);
    else
     mu = 0;
    end
    %mu=num2str(mu);
end
