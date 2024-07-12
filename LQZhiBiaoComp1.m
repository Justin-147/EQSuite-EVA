%�й�����ֵ���Ԥ���о����������ƣ�������ʱ��2017-6-2��liuqi@cea-ies.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����Ϊ������鿴ָ��Ч��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZhiBiao=LQZhiBiaoComp1(BestRR,BestThreshold1,BestThreshold2,BestHitsrate,BestR0,BestGains,BestSig,BestEnum,DisRange,DTwindow,MM,FnameQZ,pflag)
[aBestRR,bBestRR,cBestRR]=size(BestRR);
ZhiBiao=NaN(bBestRR,13);%���ָ����Ӧ����
%��Ӧʵ�ʵ�����С��2����������ǵ�����ȶ��ԣ����迼��
ind=BestEnum<2;
BestRR(ind)=NaN;
BestThreshold1(ind)=NaN;
BestThreshold2(ind)=NaN;
BestHitsrate(ind)=NaN;
BestR0(ind)=NaN;
BestGains(ind)=NaN;
BestSig(ind)=NaN;

for kk=1:1:bBestRR%����ѭ��
    tmpRR=BestRR(:,kk,:);
    if prod(isnan(tmpRR),'all')%Rֵ��Ϊ��������
        continue;
    end
    tmpRR=reshape(tmpRR,aBestRR,cBestRR);
    tmpEnum=BestEnum(:,kk,1);
    [~,innn]=unique(tmpEnum);
    innn=setdiff(innn,find(isnan(tmpEnum)));
    tmpRR=tmpRR(innn,:);%������С�ռ�  
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
    maxRR=nanmax(ZRR);%���Rֵ
    InRR=find(ZRR==maxRR);%��Ӧλ�ã����ܶ���ʱ��Ԥ����Χ��Ӧ��Rֵ��ͬ
    InRRall=InRR;
    tmpPara=NaN(length(InRR),13);
    tmpPara(:,1)=MM(kk,1);%1������
    tmpPara(:,2)=MM(kk,2);%2������
    tmpPara(:,3)=SelectDataT(BestThreshold1(:,kk,:),aBestRR,cBestRR,InRR,innn);%3������ֵ
    tmpPara(:,4)=SelectDataT(BestThreshold2(:,kk,:),aBestRR,cBestRR,InRR,innn);%4������ֵ
    tmpPara(:,5)=XDTwindow(InRR);%5Ԥ��ʱ��
    tmpPara(:,6)=YDisRange1(InRR);%6Ԥ���ռ䴰��ʼ
    tmpPara(:,7)=YDisRange(InRR);%7Ԥ���ռ䴰��ֹ
    tmpPara(:,8)=maxRR;%8Rֵ
    tmpPara(:,9)=SelectDataT(BestEnum(:,kk,:),aBestRR,cBestRR,InRR,innn);%9ʵ�ʵ�����
    tmpPara(:,10)=SelectDataT(BestHitsrate(:,kk,:),aBestRR,cBestRR,InRR,innn);%10��׼��
    tmpPara(:,11)=SelectDataT(BestGains(:,kk,:),aBestRR,cBestRR,InRR,innn);%11��������
    tmpPara(:,12)=SelectDataT(BestSig(:,kk,:),aBestRR,cBestRR,InRR,innn);%12����ˮƽ
    tmpPara(:,13)=SelectDataT(BestR0(:,kk,:),aBestRR,cBestRR,InRR,innn);%13R0
    if length(InRR)>1
        maxGain=nanmax(tmpPara(:,11));%����������
        InRR=find(tmpPara(:,11)==maxGain);%��Ӧλ��
        tmpPara=tmpPara(InRR,:);
        if length(InRR)>1
            minSig=nanmin(tmpPara(:,12));%��С����ˮƽ
            InRR=find(tmpPara(:,12)==minSig);%��Ӧλ��
            tmpPara=tmpPara(InRR,:);
            if length(InRR)>1
                minR0=nanmin(tmpPara(:,13));%��СR0����ʵ���ñ��ˣ���С����ˮƽ�����R0Ҳ���
                InRR=find(tmpPara(:,13)==minR0);%��Ӧλ��
                tmpPara=tmpPara(InRR,:);
            end
        end
    end
    ZhiBiao(kk,:)=tmpPara(end,:);%�����δ���֣�ֱ��ȡ���1������ʱ��ռ����С��
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if pflag==1
        %��ͼ����,�Ȼ��ײ�,�󸲸Ǳ��С����
        NoZS=numel(tmpRR);
        tmp1=zeros(NoZS,1);
        tmp2=flipud(reshape(repmat(DTwindow,size(tmpRR,1),1),NoZS,1));
        xx=[tmp1,tmp2,tmp2,tmp1];
        tmp1=flipud(reshape(repmat(tmpDisRange(:,1),1,size(tmpRR,2)),NoZS,1));
        tmp2=flipud(reshape(repmat(tmpDisRange(:,2),1,size(tmpRR,2)),NoZS,1));
        yy=[tmp1,tmp1,tmp2,tmp2];
        
        FNa='Times New Roman';
        FSa=10;
        FNL='��Բ';
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
        ylabel(colorbar,'Rֵ','FontName',FNL,'FontSize',FSL);
        xlabel('Ԥ��ʱ��/��','FontName',FNL,'FontSize',FSL);
        ylabel('Ԥ����Χ/km','FontName',FNL,'FontSize',FSL);
        title({['����',num2str(ZhiBiao(kk,1)),'�����ϵ���, ���������ֵ:',num2str(ZhiBiao(kk,3)),', ���������ֵ:',num2str(ZhiBiao(kk,4)),', Ԥ��ʱ��:',num2str(ZhiBiao(kk,5)),'��, Ԥ����Χ:',num2str(ZhiBiao(kk,6)),'-',num2str(ZhiBiao(kk,7)),'km'];...
            ['Rֵ:',num2str(ZhiBiao(kk,8),'%.4f'),', R0ֵ:',num2str(ZhiBiao(kk,13),'%.4f'),', ʵ�ʷ���������:',num2str(ZhiBiao(kk,9)),', ��׼��:',num2str(ZhiBiao(kk,10),'%.4f'),', ��������:',num2str(ZhiBiao(kk,11),'%.4f'),', ����ˮƽ:',num2str(ZhiBiao(kk,12))]},...
            'FontName',FNL,'FontSize',FSL);
        disp({['����',num2str(ZhiBiao(kk,1)),'�����ϵ���, ���������ֵ:',num2str(ZhiBiao(kk,3)),', ���������ֵ:',num2str(ZhiBiao(kk,4)),', Ԥ��ʱ��:',num2str(ZhiBiao(kk,5)),'��, Ԥ����Χ:',num2str(ZhiBiao(kk,6)),'-',num2str(ZhiBiao(kk,7)),'km'];...
            ['Rֵ:',num2str(ZhiBiao(kk,8),'%.4f'),', R0ֵ:',num2str(ZhiBiao(kk,13),'%.4f'),', ʵ�ʷ���������:',num2str(ZhiBiao(kk,9)),', ��׼��:',num2str(ZhiBiao(kk,10),'%.4f'),', ��������:',num2str(ZhiBiao(kk,11),'%.4f'),', ����ˮƽ:',num2str(ZhiBiao(kk,12))]});
        Figname=[FnameQZ,'_MCI',num2str(ZhiBiao(kk,1))];%�ۺ�ָ��ͼ��ÿ������ÿ���𼶵�һ��ͼ
        print(hp,[Figname,'.png'],'-dpng','-r600');
        saveas(hp,Figname,'pdf');
		close all;
    end    
end
end


function out=SelectDataT(tmpPara,aBestRR,cBestRR,InRR,innn)
tmpPara=reshape(tmpPara,aBestRR,cBestRR);
tmpPara=tmpPara(innn,:);%������С�ռ�  
tmpPara=rot90(tmpPara,2);
tmpPara=tmpPara(:);
out=tmpPara(InRR);
end


