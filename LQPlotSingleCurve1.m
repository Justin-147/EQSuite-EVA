%中国地震局地震预测研究所刘琦编制，最后调试时间2022-3-5，liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%绘制单个曲线检验结果
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutDZML=LQPlotSingleCurve1(ZhiBiao,YL,SelectEC,DisRange,Utime,Udata,FnameQZ)
Utime=Utime(:);
Udata=Udata(:);
for mm=1:1:size(ZhiBiao,1)%震级
    if isnan(ZhiBiao(mm,1))
        continue;
    end
    ind=DisRange(:,2)==ZhiBiao(mm,7);
    DZML=SelectEC(ind,mm);
    FNa='Times New Roman';
    FSa=10;
    FNL='幼圆';
    FSL=9;
    tmp1=(max(Udata)-min(Udata))/3;
    yrange=tmp1*4;
    hp=figure; hold on;
    set(hp,'Position',[507 299 630 350]);
    set(hp,'PaperPositionMode','auto');
    ob1=plot(Utime,Udata,'k');%待评估曲线
    ob2=plot([min(Utime),max(Utime)],[ZhiBiao(mm,3),ZhiBiao(mm,3)],'m:','linewidth',2);%负方向阈值
    ob3=plot([min(Utime),max(Utime)],[ZhiBiao(mm,4),ZhiBiao(mm,4)],'m:','linewidth',2);%正方向阈值
    
    ind1=Udata<ZhiBiao(mm,3)&~isnan(Udata);
    CoUdata1=Udata(ind1);
    CoUtime1=Utime(ind1);
    BaseL1=ZhiBiao(mm,3)*ones(length(CoUdata1),1);
    
    ind2=Udata>ZhiBiao(mm,4)&~isnan(Udata);
    CoUdata2=Udata(ind2);
    CoUtime2=Utime(ind2);
    BaseL2=ZhiBiao(mm,4)*ones(length(CoUdata2),1);
    
    
    ob4=fill([CoUtime1;flipud(CoUtime1)],[CoUdata1;BaseL1],'r');%负方向超阈值填充
    ob5=fill([CoUtime2;flipud(CoUtime2)],[CoUdata2;BaseL2],'r');%正方向超阈值填充
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
    tmp2=DZML.DateNum(SeDZind);%报准地震
    DAllCoUtime=diff(AllCoUtime);
    cankaoJG=DAllCoUtime(1);%参考间隔
    indJD=find(DAllCoUtime~=cankaoJG);
    if isempty(indJD)
        ShiChuang=[AllCoUtime(1);AllCoUtime(end);AllCoUtime(end);AllCoUtime(1)];
    else
        tmp3=[[AllCoUtime(1);AllCoUtime(indJD+1)],[AllCoUtime(indJD);AllCoUtime(end)]];
        ShiChuang=[tmp3,fliplr(tmp3)]';
    end
    tmp4=[(min(Udata)+tmp1*3.2);(min(Udata)+tmp1*3.2);(min(Udata)+tmp1*3.8);(min(Udata)+tmp1*3.8)];
    tmp4=repmat(tmp4,1,size(ShiChuang,2));
    ob6=fill(ShiChuang,tmp4,'m','FaceAlpha',0.4,'linewidth',0.1);%预报时窗
    ob7=scatter(DZML.DateNum,(min(Udata)+tmp1*3.5)*ones(size(DZML.DateNum)),(DZML.Magnitude-4)*50,'c','filled','markerfacecolor','b','markeredgecolor','y','linewidth',0.1);%全部地震
    ob8=scatter(tmp2,(min(Udata)+yrange*3.5/4)*ones(size(tmp2)),(DZML.Magnitude(SeDZind)-4)*50,'c','filled','markerfacecolor','r','markeredgecolor','y','linewidth',0.1);%报准地震
    if isempty(ob4)
        obb=ob5(1);
    else
        obb=ob4(1);
    end
    legend([ob1,ob2,obb,ob6(1),ob7,ob8],[{'数据曲线','异常阈值','异常时段','预报时段','漏报地震','报准地震'}],'location','eastoutside','FontName',FNL,'FontSize',FSa);
    set(gca,'position',[0.1683 0.1119 0.6349 0.7826],'FontName',FNa,'FontSize',FSa);
    xlabel('观测时间/年-月','FontName',FNL,'FontSize',FSL);
    ylabel(YL,'FontName',FNL,'FontSize',FSL);
    title({['对于',num2str(ZhiBiao(mm,1)),'级以上地震, 负向最佳阈值:',num2str(ZhiBiao(mm,3)),', 正向最佳阈值:',num2str(ZhiBiao(mm,4)),', 预报时窗:',num2str(ZhiBiao(mm,5)),'天, 预报范围:',num2str(ZhiBiao(mm,6)),'-',num2str(ZhiBiao(mm,7)),'km'];...
        ['R值:',num2str(ZhiBiao(mm,8),'%.4f'),', R0值:',num2str(ZhiBiao(mm,13),'%.4f'),', 实际发生地震数:',num2str(ZhiBiao(mm,9)),', 报准率:',num2str(ZhiBiao(mm,10),'%.4f'),', 概率增益:',num2str(ZhiBiao(mm,11),'%.4f'),', 显著水平:',num2str(ZhiBiao(mm,12))]},...
        'FontName',FNL,'FontSize',FSL);
    hold off;
    Figname=[FnameQZ,'_MSP',num2str(ZhiBiao(mm,1))];%单个曲线，每个测项每个震级档一张图
    print(hp,[Figname,'.png'],'-dpng','-r600');
    saveas(hp,Figname,'pdf');
    close all;		
end
end