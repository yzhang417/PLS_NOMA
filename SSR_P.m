clc;
clear;
close all;

%circulation
LOOP = 1000;
%Power of Noise
sigmadBm = -70;
sigma = 10^(sigmadBm/10)*10^(-3); 
%path loss factor
a = 3.0;
%dm
dm = 80;
%de
de = 80;
%M
M_array = [2 3 4];
%Q
QoS = 1.0; 
%Ptot
PtotmindBm = 0;
PtotmaxdBm = 40;
dPdBm = 2.5;

%Result Array
len = length(PtotmindBm:dPdBm:PtotmaxdBm);
Rs = zeros(LOOP,len);
AverageRs = zeros(4,len);

for index=1:1:3
    index = index
    % Number of licit user
    M = M_array(index);
    %Qm
    Q = QoS*ones(1,M);
    for loop = 1:1:LOOP
        %% Channel Information
        % Channel gain of licit user 
        h = 1/sqrt(2)*randn(1,M)+1/sqrt(2)*randn(1,M)*1j;
        h = h.*conj(h);
        nh = sort(h)./(1+dm.^a);                                    
        % Channel gain of illicit user
        he = 1/sqrt(2)*randn(1,1)+1/sqrt(2)*randn(1,1)*1j;
        nhe = he.*conj(he)/(1+de^a);      
        P = zeros(1,M);
        B = (2.^(Q)-1)./nh;
        for m = M : -1 : 1
            P(m) = B(m)*(nh(m)*sum(P(m+1:M)) + sigma); 
        end
        Pmin = sum(P);
        PmindBm = 10*log10(Pmin*10^3);  
        compter=0;
        while PmindBm > PtotmindBm
            compter = compter+1;
            if compter>20 
                index = index;
                loop = loop;
                break;
            end          
            index;
            loop = loop;
            % Channel gain of licit user 
            h = 1/sqrt(2)*randn(1,M)+1/sqrt(2)*randn(1,M)*1j;
            h = h.*conj(h);
            nh = sort(h)./(1+dm.^a);                                    
            % Channel gain of illicit user
            he = 1/sqrt(2)*randn(1,1)+1/sqrt(2)*randn(1,1)*1j;
            nhe = he.*conj(he)/(1+de^a);      
            P = zeros(1,M);
            B = (2.^(Q)-1)./nh;
            for m = M : -1 : 1
                P(m) = B(m)*(nh(m)*sum(P(m+1:M)) + sigma); 
            end
            Pmin = sum(P);
            PmindBm = 10*log10(Pmin*10^3);
        end
        Ptotindex = 0;
        for PtotdBm = PtotmindBm : dPdBm : PtotmaxdBm
            Ptotindex = Ptotindex+1;
            Ptot = 10.^(PtotdBm/10)*10^(-3);  
            ropt = zeros(1,M);
            A = (2.^(Q)-1)./(nh*Ptot);
            for m = 1 : 1 : M-1
                ropt(m) = A(m)*(Ptot*nh(m)*(1-sum(ropt(1:m-1))) + sigma)/(2^(QoS)); 
            end
            ropt(M) = 1-sum(ropt(1:M-1));

            Me = length( find( nh <= nhe ));
            Rbm = zeros(1,M);
            Rem = zeros(1,M);
            for m = Me+1:1:M
                Rbm(m) = log2( 1 + (Ptot*nh(m)*ropt(m))/(Ptot*nh(m)*sum(ropt(m+1:M))+sigma) );
                Rem(m) = log2( 1 + (Ptot*nhe*ropt(m))/(Ptot*nhe*sum(ropt(m+1:M))+sigma) );
            end
            Rsm = Rbm - Rem;
            Rs(loop,Ptotindex) = sum(Rsm(Me+1:M));
        end       
    end
    AverageRs(index,:) = nanmean(Rs);
end

for loop=1:1:LOOP
    Ptotindex = 0;
    % Channel gain of licit user 
    h = 1/sqrt(2)*randn(1,1)+1/sqrt(2)*randn(1,1)*1j;
    nh = h.*conj(h)./(1+dm.^a);                                    
    % Channel gain of illicit user
    he = 1/sqrt(2)*randn(1,1)+1/sqrt(2)*randn(1,1)*1j;
    nhe = he.*conj(he)/(1+de^a);  
    
    for PtotdBm = PtotmindBm : dPdBm : PtotmaxdBm
        Ptotindex = Ptotindex + 1;
        Ptot = 10.^(PtotdBm/10)*10^(-3);  
        Rs(loop,Ptotindex) = max(log2(1+Ptot*nh/sigma)-log2(1+Ptot*nhe/sigma),0);
    end
end
AverageRs(4,:) = nanmean(Rs);

figure;
Ptotrange = PtotmindBm:dPdBm:PtotmaxdBm;
plot(Ptotrange,AverageRs(1,:),'-*r'); hold on;
plot(Ptotrange,AverageRs(2,:),'-ob'); hold on;
plot(Ptotrange,AverageRs(3,:),'-<k'); hold on;
plot(Ptotrange,AverageRs(4,:),'-dg'); hold on;
grid on;
xlabel('P (dBm)');ylabel('Average SSR (bits/s/Hz)');
legend(['NOMA M=',num2str(M_array(1))],['NOMA M=',num2str(M_array(2))],['NOMA M=',num2str(M_array(3))], 'Conventional OMA');

