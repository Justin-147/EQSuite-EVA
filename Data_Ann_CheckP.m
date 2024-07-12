%��������Լ��
%2022-02-25������
function Ratio_Ann=Data_Ann_CheckP(data_p,ann_cyc1,FlagPanC,PlotLanMode,FF,PP)
%data_pΪ��ֵ
dt=24*60*60;%��λs
NFFT1 = 2^nextpow2(length(data_p));
fft1 = fft(data_p,NFFT1)/length(data_p);
f1 = 0.5/dt*linspace(0,1,NFFT1/2+1);
% Plot single-sided amplitude spectrum.
ann_fre=1./(ann_cyc1*24*60*60);
ind_ann=[find(f1<=ann_fre(1),1,'last'):1:find(f1>=ann_fre(2),1,'first')];
ymax_ann=2*max(abs(fft1(ind_ann)));%�������ź�������
x_ann=[ann_fre(1),ann_fre(1),ann_fre(2),ann_fre(2),ann_fre(1)];
y_ann=[0,ymax_ann,ymax_ann,0,0];
ymax_other=2*max(abs(fft1(ind_ann(end)+1:NFFT1/2+1)));%�����������ź������ȣ���������Ƶ���ƣ�
Ratio_Ann=ymax_ann/ymax_other;
if FlagPanC==1%ȫƵ�λ���
    figure;
    semilogx(f1,2*abs(fft1(1:NFFT1/2+1)),'k');
    hold on;
    semilogx(x_ann,y_ann,'r');
    hold off;
elseif FlagPanC==2%�ض�Ƶ�λ���
    figure;
    semilogx(f1(ind_ann(1):end),2*abs(fft1(ind_ann(1):NFFT1/2+1)),'k');
    hold on;
    semilogx(x_ann,y_ann,'r');
    hold off;
else
    return;
end
if PlotLanMode==1%���ı�ע
    title('FFT�׷���');
    xlabel('Ƶ�� / Hz');
    ylabel('��ֵ');
else%Ӣ�ı�ע
    title('FFT Spectrum Analysis');
    xlabel('Frequency / Hz');
    ylabel('Amplitude');
end
set(gca,'tickdir','out');
axis tight;
box off;
f_nn=find(FF=='.',1,'last')-1;
Figname=[PP,FF(1:f_nn),'_��������ȼ��'];
print(gcf,[Figname,'.png'],'-dpng','-r600');
saveas(gcf,Figname,'pdf');
close all;
end