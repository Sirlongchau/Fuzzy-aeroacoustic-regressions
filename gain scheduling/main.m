%-------------------------------------------------------------------------%
% This code is written by Hugo Henry to serve as a fuzzy gain scheduling  %
% For aero accoustic regression benchmarking for the 12/31/2025           %
%                                                                         %
% This code is aimed at research and no business usages                   %
%                                                                         %
%-------------------------------------------------------------------------%
clear; close all; clc;
%% data aquisition
% we use a regression dataset for airfoil self-noise
data_proc;
global datatrain
global membership
%% process generation
% for this research we aim to construct a fuzzy gain scheduling system to 
% approximate the Amiet analytical method through uncertainty handling 
% via a genetic fuzzy system, the Amiet model for trailing edge noise 
% generation is under the form S_pp(f)=∫Φ_pp(k_x,f)(|H(f,θ)|^2)dk_x
% we will use the symplified logarithmic version 
% SPL(f)=C+5logU+2logδ+logF(fδ/U​,α,Re) with logF≈2log(fδ/U​)−βfδ/U​+G(α,Re)
% and under a fuzzy system we will use the genetic algorithm to create the 
% constants and deal with measurement noise and uncertainty with the following form
% model becomes: SPL(f)=c_0+5logU+2logδ+c_1log(fδ/U)+c_2fδ/U+G(α,Re)
% with constants being c_0,c_1, c_2 and G

%% initialization of the training and testing datasets
train=randperm(length(processed),ceil(1503*0.8));
points=[];
for p=1:length(processed)
    if isempty(find(train==p))
        points=[points p];
    end
end

datatrain=processed(train,:);
datatest=processed(points,:);

%% initialization of the clusters
[clusters, membership,fobj] = fcm(datatrain(:,[3 5 6]), 15);
membership=membership';
%% initialisation of the GA population
popsize=150;
pop=init(popsize);

%% running optimization

options = optimoptions('ga', ...
    'Display','iter', ...
    'PlotFcn', {@gaplotbestf, @gaplotbestindiv});
chr=ga(@fitness1st,45,options);
%chr=ga(@fitness,60,options);




%% visualization 0 order clustered
% for i = 1:length(datatrain)
% 
%     c0=membership(i,1)*chr(1)+membership(i,2)*chr(2)+membership(i,3)*chr(3)+membership(i,4)*chr(4)+membership(i,5)*chr(5)+membership(i,6)*chr(6)+membership(i,7)*chr(7)+membership(i,8)*chr(8)+membership(i,9)*chr(9)+membership(i,10)*chr(10)+membership(i,11)*chr(11)+membership(i,12)*chr(12)+membership(i,13)*chr(13)+membership(i,14)*chr(14)+membership(i,15)*chr(15);
%     c1=membership(i,1)*chr(16)+membership(i,2)*chr(17)+membership(i,3)*chr(18)+membership(i,4)*chr(19)+membership(i,5)*chr(20)+membership(i,6)*chr(21)+membership(i,7)*chr(22)+membership(i,8)*chr(23)+membership(i,9)*chr(24)+membership(i,10)*chr(25)+membership(i,11)*chr(26)+membership(i,12)*chr(27)+membership(i,13)*chr(28)+membership(i,14)*chr(29)+membership(i,15)*chr(30);
%     c2=membership(i,1)*chr(31)+membership(i,2)*chr(32)+membership(i,3)*chr(33)+membership(i,4)*chr(34)+membership(i,5)*chr(35)+membership(i,6)*chr(36)+membership(i,7)*chr(37)+membership(i,8)*chr(38)+membership(i,9)*chr(39)+membership(i,10)*chr(40)+membership(i,11)*chr(41)+membership(i,12)*chr(42)+membership(i,13)*chr(43)+membership(i,14)*chr(44)+membership(i,15)*chr(45);
%     G=membership(i,1)*chr(46)+membership(i,2)*chr(47)+membership(i,3)*chr(48)+membership(i,4)*chr(49)+membership(i,5)*chr(50)+membership(i,6)*chr(51)+membership(i,7)*chr(52)+membership(i,8)*chr(53)+membership(i,9)*chr(54)+membership(i,10)*chr(55)+membership(i,11)*chr(56)+membership(i,12)*chr(57)+membership(i,13)*chr(58)+membership(i,14)*chr(59)+membership(i,15)*chr(60);
% 
%     outputtrain(i)=c0+5*log10(datatrain(i,4))+2*log10(datatrain(i,6))+c1*log10(datatrain(i,1)*datatrain(i,6)/datatrain(i,4))+c2*datatrain(i,1)*datatrain(i,6)/datatrain(i,4)+G ;
% end
% fittrain=sum(outputtrain-processed(:,end));
% 
% for i = 1:length(datatest)
%     c0=membership(i,1)*chr(1)+membership(i,2)*chr(2)+membership(i,3)*chr(3)+membership(i,4)*chr(4)+membership(i,5)*chr(5)+membership(i,6)*chr(6)+membership(i,7)*chr(7)+membership(i,8)*chr(8)+membership(i,9)*chr(9)+membership(i,10)*chr(10)+membership(i,11)*chr(11)+membership(i,12)*chr(12)+membership(i,13)*chr(13)+membership(i,14)*chr(14)+membership(i,15)*chr(15);
%     c1=membership(i,1)*chr(16)+membership(i,2)*chr(17)+membership(i,3)*chr(18)+membership(i,4)*chr(19)+membership(i,5)*chr(20)+membership(i,6)*chr(21)+membership(i,7)*chr(22)+membership(i,8)*chr(23)+membership(i,9)*chr(24)+membership(i,10)*chr(25)+membership(i,11)*chr(26)+membership(i,12)*chr(27)+membership(i,13)*chr(28)+membership(i,14)*chr(29)+membership(i,15)*chr(30);
%     c2=membership(i,1)*chr(31)+membership(i,2)*chr(32)+membership(i,3)*chr(33)+membership(i,4)*chr(34)+membership(i,5)*chr(35)+membership(i,6)*chr(36)+membership(i,7)*chr(37)+membership(i,8)*chr(38)+membership(i,9)*chr(39)+membership(i,10)*chr(40)+membership(i,11)*chr(41)+membership(i,12)*chr(42)+membership(i,13)*chr(43)+membership(i,14)*chr(44)+membership(i,15)*chr(45);
%     G=membership(i,1)*chr(46)+membership(i,2)*chr(47)+membership(i,3)*chr(48)+membership(i,4)*chr(49)+membership(i,5)*chr(50)+membership(i,6)*chr(51)+membership(i,7)*chr(52)+membership(i,8)*chr(53)+membership(i,9)*chr(54)+membership(i,10)*chr(55)+membership(i,11)*chr(56)+membership(i,12)*chr(57)+membership(i,13)*chr(58)+membership(i,14)*chr(59)+membership(i,15)*chr(60);
% 
%     outputtest(i)=c0+5*log10(datatest(i,4))+2*log10(datatest(i,6))+c1*log10(datatest(i,1)*datatest(i,6)/datatest(i,4))+c2*datatest(i,1)*datatest(i,6)/datatest(i,4)+G ;
% end
% fittest=sum(outputtest-processed(:,end));
%% visu first order tsk
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

        outputtrain(i,1) = ...
            5*log10(datatrain(i,4)) + ...
            2*log10(datatrain(i,6)) + ...
            y_fuzzy;
end
for i=1:length(datatest)
        z1 = log10(datatest(i,1)*datatest(i,6)/datatest(i,4));
        z2 = (datatest(i,1)*datatest(i,6)/datatest(i,4));

        y_fuzzy = 0;

            for k = 1:15
                a0 = chr(3*(k-1)+1);
                a1 = chr(3*(k-1)+2);
                a2 = chr(3*(k-1)+3);

                y_fuzzy = y_fuzzy + membership(i,k) * ...
                    (a0 + a1*z1 + a2*z2);
            end

        outputtest(i,1) = ...
            5*log10(datatest(i,4)) + ...
            2*log10(datatest(i,6)) + ...
            y_fuzzy;
end


figure;
title("Training set")
scatter(train,outputtrain,"red")
hold on
scatter(train,datatrain(:,end),"blue")
hold off

figure;
title("Testing set")
scatter(points,outputtest,"red")
hold on
scatter(points,datatest(:,end),"blue")
hold off