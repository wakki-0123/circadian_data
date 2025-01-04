%fs=250;
fs=250;
%fs=70;
epoch_l=20;
ch_size=8;

e_tail=2500;

% For MSE

m=2;
r=0.2;
factor=100;
size_sur_data=10;
maxiter=100;

ch_nmae={'F3','F4','C3','C4','P3','P4','O1','O2'};

normal_sub_list=importdata('sub_list_normal.txt');
asp_sub_list=importdata('sub_list_asp.txt');

addpath('./EntropyHub_v2.0.0/');



for s_i=1:1:size(normal_sub_list,1)
    s_i
    load_data=importdata(char(normal_sub_list(s_i)),',',1);
    time_series=load_data.data;


    
    for ch=1:1:ch_size
        ch
        passband = [0.1/(fs/2) 30/(fs/2)];
        
        fir = fir1 ( floor ( size(time_series(1:end-e_tail,ch+1),1) / 3 ) - 1, passband );
        
        filtered_ch = filtfilt(fir, 1, double(time_series(1:end-e_tail,ch+1)));

        
        for epoch_i=1:1:floor(length(filtered_ch)/(epoch_l*fs))
            epoch_i
            %[psdx(s_i,ch,epoch_i,:),freq]=pwelch((time_series(:,ch+1)),(epoch_l*fs),(epoch_l*fs)/2,[0:0.5:40],fs);
            e_s=(epoch_i-1)*fs*epoch_l+1;
            e_e=(epoch_i)*fs*epoch_l;

            Mobj = MSobject('FuzzEn');
            MSx=cMSEn(filtered_ch(e_s:e_e),Mobj,'Scales',factor);
            
            HC_MSE(s_i,ch,epoch_i,:)=MSx;
            
        end
    end
        
end


for s_i=1:1:size(asp_sub_list,1)
    s_i
    load_data=importdata(char(asp_sub_list(s_i)),',',1);
    time_series=load_data.data;


    
    for ch=1:1:ch_size
        ch
        passband = [0.1/(fs/2) 30/(fs/2)];        
        fir = fir1 ( floor ( size(time_series(1:end-e_tail,ch+1),1) / 3 ) - 1, passband);
        filtered_ch = filtfilt(fir, 1, double(time_series(1:end-e_tail,ch+1)));

        

        for epoch_i=1:1:floor(length(filtered_ch)/(epoch_l*fs))
            epoch_i
            %[psdx(s_i,ch,epoch_i,:),freq]=pwelch((time_series(:,ch+1)),(epoch_l*fs),(epoch_l*fs)/2,[0:0.5:40],fs);
            e_s=(epoch_i-1)*fs*epoch_l+1;
            e_e=(epoch_i)*fs*epoch_l;

            Mobj = MSobject('FuzzEn');
            MSx=cMSEn(filtered_ch(e_s:e_e),Mobj,'Scales',factor);

            AS_MSE(s_i,ch,epoch_i,:)=MSx;

            
        end
    end
        
end



save('MSE_result_fuzzy.mat');


for ch=1:1:ch_size

    health_m=reshape(mean(mean(HC_MSE(:,ch,:,:),3),1),[1 factor]);
    health_ste=reshape(std(mean(HC_MSE(:,ch,:,:),3))/sqrt(size(normal_sub_list,1)),[1 factor]);

    iaaft_health_m=reshape(mean(mean(iaaft_HC_MSE(:,ch,:,:),3),1),[1 factor]);
    iaaft_health_ste=reshape(std(mean(iaaft_HC_MSE(:,ch,:,:),3))/sqrt(size(normal_sub_list,1)),[1 factor]);

    



    
    figure;

    plot([1:1:factor],health_m,'-k','Linewidth',2);

    hold on;
    
    plot([1:1:factor],health_m+health_ste,':k','Linewidth',2);
    plot([1:1:factor],health_m-health_ste,':k','Linewidth',2);


    plot([1:1:factor],iaaft_health_m,'-r','Linewidth',2);
    
    plot([1:1:factor],iaaft_health_m+iaaft_health_ste,':r','Linewidth',2);
    plot([1:1:factor],iaaft_health_m-iaaft_health_ste,':r','Linewidth',2);


    plot([1:1:factor],reshape(mean(AS_MSE(1,ch,:,:),3),[1 factor]),'-b','Linewidth',4);
    plot([1:1:factor],reshape(mean(AS_MSE(2,ch,:,:),3),[1 factor]),'-c','Linewidth',4);

    
    %semilogy(freq,reshape(mean(as_psdx(1,ch,:,:),3),[1 81]),'-r','Linewidth',2);
    %semilogy(freq,reshape(mean(as_psdx(2,ch,:,:),3),[1 81]),'-m','Linewidth',2);
    
    grid on;
    title(char(ch_nmae(ch)));


    xlabel('time-scale')
    ylabel('SampEn');
    
    
    saveas(gcf,strcat('mse_',num2str(ch)),'epsc');
end
