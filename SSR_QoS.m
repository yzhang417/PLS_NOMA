clc;
close all;

Index=1000;
d=80;
a=3;
sigma_dbm=-70;
sigma=10^((sigma_dbm-30)/10);
M_array=[2,3,4];
M_len=length(M_array);
P_dbm=20;
P=10^((P_dbm-30)/10);
Q_min=0;
Q_max=6;
d_Q=0.25;
Q_range=Q_min:d_Q:Q_max;
len=length(Q_range);
result_array=zeros(Index,len,4);
average_result_array=zeros(4,len);

for M_index=1:1:M_len
    M=M_array(M_index);
    for index=1:1:Index
        gain=sqrt(1/2)*randn(1,M)+sqrt(1/2)*randn(1,M)*1i;
        h=gain./d^(a/2);
        h1=h.*conj(h);
        h2=sort(h1);
        gain_eaves=sqrt(1/2)*randn+sqrt(1/2)*randn*1i;
        h_eaves=gain_eaves/d^(a/2);
        h_eaves=h_eaves*conj(h_eaves);
        for Q_index=1:1:len
            Q=Q_range(Q_index);
            B=(2.^(Q)-1)./h2;
            P_feasibility=zeros(1,M);
            for m=M:-1:1        
                P_feasibility(m)=B(m)*(h2(m)*sum(P_feasibility(m+1:M))+sigma); 
            end
            P_feas=sum(P_feasibility);
            P_feas_dbm=10*log10(P_feas)+30;
            if P_dbm<P_feas_dbm
                result_array(index,Q_index,M_index)=0;
            else
                A=(2.^(Q)-1)./(P*h2);
                r=zeros(1,M);
                for m=1:1:M-1
                    r(m)=(A(m)*(P*h2(m)*(1-sum(r(1:m-1)))+sigma))./(2.^Q);
                end
                r(M)=1-sum(r(1:M-1));
                R_user=zeros(1,M);
                R_eaves=zeros(1,M);
                R_secrecy=zeros(1,M);
                for m=1:1:M
                    R_user(m)=log2(1+P*h2(m)*r(m)/(P*h2(m)*sum(r(m+1:M))+sigma));
                    R_eaves(m)=log2(1+P*h_eaves*r(m)/(P*h_eaves*sum(r(m+1:M))+sigma));
                    if R_user(m)<R_eaves(m)
                        R_secrecy(m)=0;
                    else
                        R_secrecy(m)=R_user(m)-R_eaves(m);
                    end
                end
                R_sec=sum(R_secrecy);
                result_array(index,Q_index,M_index)=R_sec;
            end
        end
    end
    average_result_array(M_index,:)=nanmean(result_array(:,:,M_index));
end

for index=1:1:Index
    gain=sqrt(1/2)*randn+sqrt(1/2)*randn*1i;
    h=gain/d^(a/2);
    h1=h*conj(h); 
    gain_eaves=sqrt(1/2)*randn+sqrt(1/2)*randn*1i;
    h_eaves=gain_eaves/d^(a/2);
    h_eaves=h_eaves*conj(h_eaves);
    for Q_index=1:1:len
        R_b=log2(1+P*h1/sigma);
        R_e=log2(1+P*h_eaves/sigma);
        if R_b<R_e
            R_secrecy=0;
        else
            R_secrecy=R_b-R_e;
        end
        result_array(index,Q_index,4)=R_secrecy;
    end
end
average_result_array(4,:)=nanmean(result_array(:,:,4));


figure;
plot(Q_range,average_result_array(1,:),'-*k'); hold on;
plot(Q_range,average_result_array(2,:),'-ob'); hold on;
plot(Q_range,average_result_array(3,:),'-sr'); hold on; 
plot(Q_range,average_result_array(4,:),'-xg'); hold on;
grid on;
xlabel('Q(bits/s/Hz)');ylabel('Average SSR (bits/s/Hz)');
legend('NOMA,M=2','NOMA,M=3','NOMA,M=4','OMA(TDMA)');
