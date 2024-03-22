%  rng default
%    %daytime_num=datenum(daytime_hikiti_heartrate,'yyyy/mm/dd HH:MM:SS');
%    %xq=(datetime(2023,07,11,13,55,00):seconds(5):datetime(2023,08,01,21,33,00));
%   x=interpft(hikiti_heartrate_AllDay,368377);
 %figure;
 %plot(timedate, heartrate,'k-')
%  figure;
%  plot(daytime_num, heartrate, 'k-')
%  datetick('x','yyyy/mm/dd HH:MM:SS','keeplimits','keepticks')
%  plot(xq,x,'k-')
% figure;
% plot(daytime_num,heartrate,'k--');
% datetick('x','yyyy/mm/dd HH:MM:SS','keeplimits','keepticks')
% fs=100;
% [pxx1,f]=pwelch(heartrate,500,300,500,fs)
% [pxx2,f]=pwelch(x,500,300,500,fs)
% figure;
% loglog(f,pxx1,'k')
% hold on
% loglog(f,pxx2,'b')
% xlabel('Frequency(Hz)')
% ylabel('PSD(dB/Hz)')


Z_maru1=zscore(maruyama_heartrate_AllDay_Resample_Frequency)