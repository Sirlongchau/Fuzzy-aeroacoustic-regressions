data=load("airfoil_self_noise.dat");

%data have : {F AoA c U Re/nu delta dB}
processed=[data(:,1) data(:,2) data(:,3) data(:,4) data(:,3).*data(:,4) data(:,5) data(:,6)]; %convert velocity and chord lentgh into normalized by nu reynolds numbers and take out the output

