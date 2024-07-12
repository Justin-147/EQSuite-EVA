function Plot_Ann(time_p,data_p,data_mean_ST_ann,data_mean_ST_other,ann_cyc1,ann_cyc2,PlotLanMode)
FSA=10;
FSL=12;
FNA='Times New Roman';
FNL='����';
nn=3;
datex=datenum(num2str(time_p),'yyyymmdd');
if PlotLanMode==1
    yla1='ԭʼ�۲�';
    yla2={'��Ƶ�ι�һ��ƽ������ANA';['(����',num2str(ann_cyc1(2)),'-',num2str(ann_cyc1(1)),'��)']};
    yla3={'����Ƶ�ι�һ��ƽ������ONA';['(����',num2str(ann_cyc2(2)),'-',num2str(ann_cyc2(1)),'��)']};
    xla='�۲����� / ��';
else
    yla1='Observation';
    yla2={'ANA';['(Period from ',num2str(ann_cyc1(2)),' to ',num2str(ann_cyc1(1)),' days)']};
    yla3={'ONA';['(Period from ',num2str(ann_cyc2(2)),' to ',num2str(ann_cyc2(1)),' days)']};
    xla='Date / Year';
end
hf=figure();
set(hf,'position',[680 202 710 776]);
set(hf,'PaperPositionMode','auto');
subplot(3*nn,1,1:nn)
plot(datex,data_p,'k');
xlim([datex(1) datex(end)])
tlabel();
set(gca,'tickdir','out');
set(gca,'xticklabel',[]);
set(gca,'FontSize',FSA,'FontName',FNA);
ylabel(yla1,'FontSize',FSL,'FontName',FNL);
box off;
subplot(3*nn,1,nn+1:2*nn)
plot(datex,data_mean_ST_ann,'k');
xlim([datex(1) datex(end)])
tlabel();
set(gca,'tickdir','out');
set(gca,'xticklabel',[]);
set(gca,'FontSize',FSA,'FontName',FNA);
ylabel(yla2,'FontSize',FSL,'FontName',FNL);
box off;
subplot(3*nn,1,2*nn+1:3*nn)
plot(datex,data_mean_ST_other,'k');
xlim([datex(1) datex(end)])
tlabel();
set(gca,'tickdir','out');
set(gca,'FontSize',FSA,'FontName',FNA);
ylabel(yla3,'FontSize',FSL,'FontName',FNL);
xlabel(xla,'FontSize',FSL,'FontName',FNL);
box off;