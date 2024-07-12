%中国地震局地震预测研究所刘琦编制，最后调试时间2017-6-2，liuqi@cea-ies.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%以震级为档计算查看指标效能
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZhiBiao=LQZhiBiaoComp1(BestRR,BestThreshold1,BestThreshold2,BestHitsrate,BestR0,BestGains,BestSig,BestEnum,DisRange,DTwindow,MM,FnameQZ,pflag)
[aBestRR,bBestRR,cBestRR]=size(BestRR);
ZhiBiao=NaN(bBestRR,13);%最佳指标相应参数
%对应实际地震数小于2的情况，考虑到结果稳定性，不予考虑
ind=BestEnum<2;
BestRR(ind)=NaN;
BestThreshold1(ind)=NaN;
BestThreshold2(ind)=NaN;
BestHitsrate(ind)=NaN;
BestR0(ind)=NaN;
BestGains(ind)=NaN;
BestSig(ind)=NaN;

for kk=1:1:bBestRR%按震级循环
    tmpRR=BestRR(:,kk,:);
    if prod(isnan(tmpRR),'all')%R值都为空则跳过
        continue;
    end
    tmpRR=reshape(tmpRR,aBestRR,cBestRR);
    tmpEnum=BestEnum(:,kk,1);
    [~,innn]=unique(tmpEnum);
    innn=setdiff(innn,find(isnan(tmpEnum)));
    tmpRR=tmpRR(innn,:);%保留最小空间  
    tmpRR=rot90(tmpRR,2);
    
    tmpDisRange=DisRange;
    tmpDisRange=tmpDisRange(innn,:); 
    YDisRange=repmat(flipud(tmpDisRange(:,2)),1,length(DTwindow));
    YDisRange=YDisRange(:);
    YDisRange1=repmat(flipud(tmpDisRange(:,1)),1,length(DTwindow));
    YDisRange1=YDisRange1(:);
    XDTwindow=repmat(fliplr(DTwindow),size(tmpDisRange,1),1);
    XDTwindow=XDTwindow(:);
    ZRR=tmpRR(:);
    maxRR=nanmax(ZRR);%最大R值
    InRR=find(ZRR==maxRR);%对应位置，可能多种时空预报范围对应的R值相同
    InRRall=InRR;
    tmpPara=NaN(length(InRR),13);
    tmpPara(:,1)=MM(kk,1);%1震级下限
    tmpPara(:,2)=MM(kk,2);%2震级上限
    tmpPara(:,3)=SelectDataT(BestThreshold1(:,kk,:),aBestRR,cBestRR,InRR,innn);%3负向阈值
    tmpPara(:,4)=SelectDataT(BestThreshold2(:,kk,:),aBestRR,cBestRR,InRR,innn);%4正向阈值
    tmpPara(:,5)=XDTwindow(InRR);%5预报时窗
    tmpPara(:,6)=YDisRange1(InRR);%6预报空间窗起始
    tmpPara(:,7)=YDisRange(InRR);%7预报空间窗终止
    tmpPara(:,8)=maxRR;%8R值
    tmpPara(:,9)=SelectDataT(BestEnum(:,kk,:),aBestRR,cBestRR,InRR,innn);%9实际地震数
    tmpPara(:,10)=SelectDataT(BestHitsrate(:,kk,:),aBestRR,cBestRR,InRR,innn);%10报准率
    tmpPara(:,11)=SelectDataT(BestGains(:,kk,:),aBestRR,cBestRR,InRR,innn);%11概率增益
    tmpPara(:,12)=SelectDataT(BestSig(:,kk,:),aBestRR,cBestRR,InRR,innn);%12显著水平
    tmpPara(:,13)=SelectDataT(BestR0(:,kk,:),aBestRR,cBestRR,InRR,innn);%13R0
    if length(InRR)>1
        maxGain=nanmax(tmpPara(:,11));%最大概率增益
        InRR=find(tmpPara(:,11)==maxGain);%对应位置
        tmpPara=tmpPara(InRR,:);
        if length(InRR)>1
            minSig=nanmin(tmpPara(:,12));%最小显著水平
            InRR=find(tmpPara(:,12)==minSig);%对应位置
            tmpPara=tmpPara(InRR,:);
            if length(InRR)>1
                minR0=nanmin(tmpPara(:,13));%最小R0，其实不用比了，最小显著水平相等则R0也相等
                InRR=find(tmpPara(:,13)==minR0);%对应位置
                tmpPara=tmpPara(InRR,:);
            end
        end
    end
    ZhiBiao(kk,:)=tmpPara(end,:);%如果还未区分，直接取最后1个，即时空占用最小的
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if pflag==1
        %绘图部分,先画底层,后覆盖表层小区域
        NoZS=numel(tmpRR);
        tmp1=zeros(NoZS,1);
        tmp2=flipud(reshape(repmat(DTwindow,size(tmpRR,1),1),NoZS,1));
        xx=[tmp1,tmp2,tmp2,tmp1];
        tmp1=flipud(reshape(repmat(tmpDisRange(:,1),1,size(tmpRR,2)),NoZS,1));
        tmp2=flipud(reshape(repmat(tmpDisRange(:,2),1,size(tmpRR,2)),NoZS,1));
        yy=[tmp1,tmp1,tmp2,tmp2];
        
        FNa='Times New Roman';
        FSa=10;
        FNL='幼圆';
        FSL=9;
        hp=figure; hold on;
        set(hp,'Position',[507 299 630 450]);
        set(hp,'PaperPositionMode','auto');
        tmp=tmpRR(:);
        fill(xx',yy',tmp');
%         indd=isnan(tmp);
%         fill(xx(indd,:)',yy(indd,:)','w');
        xlim([0 DTwindow(end)]);
        ylim([min(DisRange(:,1)) max(DisRange(:,2))]);
        colormap('jet')
        caxis([0 max(tmp)]);
        %caxis([0 1]);
        colorbar('EastOutside','FontName',FNa,'FontSize',FSa);
        hold on;
        plot(XDTwindow(InRRall),YDisRange(InRRall),'kp','markersize',9,'markerfacecolor','w');
        plot(ZhiBiao(kk,5),ZhiBiao(kk,7),'kp','markersize',12,'markerfacecolor','w');
        hold off;
        set(gca,'position',[0.1 0.10 0.75 0.79],'FontName',FNa,'FontSize',FSa);
        ylabel(colorbar,'R值','FontName',FNL,'FontSize',FSL);
        xlabel('预报时窗/天','FontName',FNL,'FontSize',FSL);
        ylabel('预报范围/km','FontName',FNL,'FontSize',FSL);
        title({['对于',num2str(ZhiBiao(kk,1)),'级以上地震, 负向最佳阈值:',num2str(ZhiBiao(kk,3)),', 正向最佳阈值:',num2str(ZhiBiao(kk,4)),', 预报时窗:',num2str(ZhiBiao(kk,5)),'天, 预报范围:',num2str(ZhiBiao(kk,6)),'-',num2str(ZhiBiao(kk,7)),'km'];...
            ['R值:',num2str(ZhiBiao(kk,8),'%.4f'),', R0值:',num2str(ZhiBiao(kk,13),'%.4f'),', 实际发生地震数:',num2str(ZhiBiao(kk,9)),', 报准率:',num2str(ZhiBiao(kk,10),'%.4f'),', 概率增益:',num2str(ZhiBiao(kk,11),'%.4f'),', 显著水平:',num2str(ZhiBiao(kk,12))]},...
            'FontName',FNL,'FontSize',FSL);
        disp({['对于',num2str(ZhiBiao(kk,1)),'级以上地震, 负向最佳阈值:',num2str(ZhiBiao(kk,3)),', 正向最佳阈值:',num2str(ZhiBiao(kk,4)),', 预报时窗:',num2str(ZhiBiao(kk,5)),'天, 预报范围:',num2str(ZhiBiao(kk,6)),'-',num2str(ZhiBiao(kk,7)),'km'];...
            ['R值:',num2str(ZhiBiao(kk,8),'%.4f'),', R0值:',num2str(ZhiBiao(kk,13),'%.4f'),', 实际发生地震数:',num2str(ZhiBiao(kk,9)),', 报准率:',num2str(ZhiBiao(kk,10),'%.4f'),', 概率增益:',num2str(ZhiBiao(kk,11),'%.4f'),', 显著水平:',num2str(ZhiBiao(kk,12))]});
        Figname=[FnameQZ,'_MCI',num2str(ZhiBiao(kk,1))];%综合指标图，每个测项每个震级档一张图
        print(hp,[Figname,'.png'],'-dpng','-r600');
        saveas(hp,Figname,'pdf');
		close all;
    end    
end
end


function out=SelectDataT(tmpPara,aBestRR,cBestRR,InRR,innn)
tmpPara=reshape(tmpPara,aBestRR,cBestRR);
tmpPara=tmpPara(innn,:);%保留最小空间  
tmpPara=rot90(tmpPara,2);
tmpPara=tmpPara(:);
out=tmpPara(InRR);
end


