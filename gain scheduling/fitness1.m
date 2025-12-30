function fit=fitness1(chr)
%chr first 24 lines are membership values and the rest are the rules (36
%rules)
membernum=5;
global datatrain
data=datatrain;
    for ind=1:length(data)
        for k=1:membernum % goes through the 5 MFs
            mu(1,k)=triangle(data(ind,1),chr(k),chr(k+1),chr(k+2));
            mu(2,k)=triangle(data(ind,1),chr(k+membernum),chr(k+membernum+1),chr(k+membernum+2));
            mu(3,k)=triangle(data(ind,1),chr(k+2*membernum),chr(k+2*membernum+1),chr(k+2*membernum+2));
            mu(4,k)=triangle(data(ind,1),chr(k+3*membernum),chr(k+3*membernum+1),chr(k+3*membernum+2));
            mu(5,k)=triangle(data(ind,1),chr(k+4*membernum),chr(k+4*membernum+1),chr(k+4*membernum+2));
        end
    
        rindex=0; % index of the rule (1 to 3125 )

        % for index=1:3125
        %     if index<=3125/5
        %         rmu(index)=
            for k=1:membernum % goes through each membership of input 1
                for k2=1:membernum
                    for k3=1:membernum
                        for k4=1:membernum
                            for k5=1:membernum% goes through each membership of input 2
                                rindex=rindex+1; % go to next rule (should go up to 9)
                                mflist=[mu(1,k) mu(2,k2) mu(3,k3) mu(4,k4) mu(5,k5)];
                                rmu(rindex)=min(mflist); % get the activation of the rule based on every combinaisons   
                            end
                        end
                    end                
                end
            end
            agregg=sum(rmu); % get the aggregation of the activation
            for p=1:rindex % go through every rule and divide by the agreggation
                output_raw(ind,p)=(rmu(p)*(chr(76+(p-1)*6)+chr(76+(p-1)*6+1)*data(ind,1)+chr(76+(p-1)*6+2)*data(ind,2) + chr(76+(p-1)*6+3)* data(ind,3) +chr(76+(p-1)*6+4)* data(ind,4) +chr(76+(p-1)*6+5)* data(ind,5)))/agregg; % i*9-9+p to get to the rule in question and +24 to bypass the 24 membership lines, only first component because TSK 0 order 
            end
        output(ind)=sum(output_raw(ind,:));
    end
    
    newout=output';
    for val=1:length(data)
        error(val)=abs(newout(val)-data(val,6));
    end
    q = 0.95;
    J3 = quantile(error, q);
    
    fit = sum(error) + 50*J3;
    f=sum(error);
end