%�й�����ֵ���Ԥ���о����������ƣ�������ʱ��2022-3-5��liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���Ƶ������߼�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutDZML=LQPlotSingleCurve1(ZhiBiao,YL,SelectEC,DisRange,Utime,Udata,FnameQZ)
Utime=Utime(:);
Udata=Udata(:);
for mm=1:1:size(ZhiBiao,1)%��
    if isnan(ZhiBiao(mm,1))
        continue;
    end
    ind=DisRange(:,2)==ZhiBiao(mm,7);
    DZML=SelectEC(ind,mm);
    FNa='Times New Roman';
    FSa=10;
    FNL='��Բ';
    FSL=9;
    tmp1=(max(Udata)-min(Udata))/3;
    yrange=tmp1*4;
    hp=figure; hold on;
    set(hp,'Position',[507 299 630 350]);
    set(hp,'PaperPositionMode','auto');
    ob1=plot(Utime,Udata,'k');%����������
    ob2=plot([min(Utime),max(Utime)],[ZhiBiao(mm,3),ZhiBiao(mm,3)],'m:','linewidth',2);%��������ֵ
    ob3=plot([min(Utime),max(Utime)],[ZhiBiao(mm,4),ZhiBiao(mm,4)],'m:','linewidth',2);%��������ֵ
    
    ind1=Udata<ZhiBiao(mm,3)&~isnan(Udata);
    CoUdata1=Udata(ind1);
    CoUtime1=Utime(ind1);
    BaseL1=ZhiBiao(mm,3)*ones(length(CoUdata1),1);
    
    ind2=Udata>ZhiBiao(mm,4)&~isnan(Udata);
    CoUdata2=Udata(ind2);
    CoUtime2=Utime(ind2);
    BaseL2=ZhiBiao(mm,4)*ones(length(CoUdata2),1);
    
    
    ob4=fill([CoUtime1;flipud(CoUtime1)],[CoUdata1;BaseL1],'r');%��������ֵ���
    ob5=fill([CoUtime2;flipud(CoUtime2)],[CoUdata2;BaseL2],'r');%��������ֵ���
    datetick('x','yyyy-mm');
    ylim([min(Udata) min(Udata)+yrange]);
    xlim([min(Utime) max(Utime)]);
    
    CoUtime=union(CoUtime1,CoUtime2);
    AllCoUtime=repmat(CoUtime,1,ZhiBiao(mm,5))+repmat([1:1:ZhiBiao(mm,5)],length(CoUtime),1);
    AllCoUtime=intersect(unique(AllCoUtime),Utime);
    SeDZind=ismember(DZML.DateNum,AllCoUtime);
    DZML.PredCorrectInd=SeDZind;
    DZML.ZhibiaoIndex=mm;
    OutDZML(mm)=DZML;
    tmp2=DZML.DateNum(SeDZind);%��׼����
    DAllCoUtime=diff(AllCoUtime);
    cankaoJG=DAllCoUtime(1);%�ο����
    indJD=find(DAllCoUtime~=cankaoJG);
    if isempty(indJD)
        ShiChuang=[AllCoUtime(1);AllCoUtime(end);AllCoUtime(end);AllCoUtime(1)];
    else
        tmp3=[[AllCoUtime(1);AllCoUtime(indJD+1)],[AllCoUtime(indJD);AllCoUtime(end)]];
        ShiChuang=[tmp3,fliplr(tmp3)]';
    end
    tmp4=[(min(Udata)+tmp1*3.2);(min(Udata)+tmp1*3.2);(min(Udata)+tmp1*3.8);(min(Udata)+tmp1*3.8)];
    tmp4=repmat(tmp4,1,size(ShiChuang,2));
    ob6=fill(ShiChuang,tmp4,'m','FaceAlpha',0.4,'linewidth',0.1);%Ԥ��ʱ��
    ob7=scatter(DZML.DateNum,(min(Udata)+tmp1*3.5)*ones(size(DZML.DateNum)),(DZML.Magnitude-4)*50,'c','filled','markerfacecolor','b','markeredgecolor','y','linewidth',0.1);%ȫ������
    ob8=scatter(tmp2,(min(Udata)+yrange*3.5/4)*ones(size(tmp2)),(DZML.Magnitude(SeDZind)-4)*50,'c','filled','markerfacecolor','r','markeredgecolor','y','linewidth',0.1);%��׼����
    if isempty(ob4)
        obb=ob5(1);
    else
        obb=ob4(1);
    end
    legend([ob1,ob2,obb,ob6(1),ob7,ob8],[{'��������','�쳣��ֵ','�쳣ʱ��','Ԥ��ʱ��','©������','��׼����'}],'location','eastoutside','FontName',FNL,'FontSize',FSa);
    set(gca,'position',[0.1683 0.1119 0.6349 0.7826],'FontName',FNa,'FontSize',FSa);
    xlabel('�۲�ʱ��/��-��','FontName',FNL,'FontSize',FSL);
    ylabel(YL,'FontName',FNL,'FontSize',FSL);
    title({['����',num2str(ZhiBiao(mm,1)),'�����ϵ���, ���������ֵ:',num2str(ZhiBiao(mm,3)),', ���������ֵ:',num2str(ZhiBiao(mm,4)),', Ԥ��ʱ��:',num2str(ZhiBiao(mm,5)),'��, Ԥ����Χ:',num2str(ZhiBiao(mm,6)),'-',num2str(ZhiBiao(mm,7)),'km'];...
        ['Rֵ:',num2str(ZhiBiao(mm,8),'%.4f'),', R0ֵ:',num2str(ZhiBiao(mm,13),'%.4f'),', ʵ�ʷ���������:',num2str(ZhiBiao(mm,9)),', ��׼��:',num2str(ZhiBiao(mm,10),'%.4f'),', ��������:',num2str(ZhiBiao(mm,11),'%.4f'),', ����ˮƽ:',num2str(ZhiBiao(mm,12))]},...
        'FontName',FNL,'FontSize',FSL);
    hold off;
    Figname=[FnameQZ,'_MSP',num2str(ZhiBiao(mm,1))];%�������ߣ�ÿ������ÿ���𼶵�һ��ͼ
    print(hp,[Figname,'.png'],'-dpng','-r600');
    saveas(hp,Figname,'pdf');
    close all;		
end
end