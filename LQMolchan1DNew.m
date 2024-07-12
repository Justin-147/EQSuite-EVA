%�й�����ֵ���Ԥ���о����������ƣ�������ʱ��2022-3-4��liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ü���Molchan���㺯��������ָ�������µĽ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LQMolchan1DNew(ZhiBiao,SelectEC,Utime,Udata,NoThreshold,DisRange,FnameQZ)
for mm=1:1:size(ZhiBiao,1)%��
    if isnan(ZhiBiao(mm,1))
        continue;
    end
    ind=DisRange(:,2)==ZhiBiao(mm,7);
    DZML=SelectEC(ind,mm);
    [~,~,N,Hits,Coverrate,AreaM]=LQMolchanCompute2_ParNew(Utime,Udata,ZhiBiao(mm,5),NoThreshold,DZML.DateNum); 
    RR=Hits/N-Coverrate;
    Dists=sqrt((1-Hits/N).^2+Coverrate.^2);
    Missrate=1-Hits/N;
    Gains=Hits/N./Coverrate;
    [~,indRR]=max(RR);
    [~,indDists]=min(Dists);
    indRR=indRR(1);
    indDists=indDists(1);
    %%%%%%%%%%%%%%%%%%%%%%
    %����
    LW1=2;
    LW2=1.2;
    LW3=2;
    FSa=10;
    FSl=9;
    FNa='times new roman';
    FNl='����';
    hn0=1:N;
    h=0:N;
    vn0=1-hn0/N;
    tao1=[];
    tao2=[];
    for ii=1:1:length(hn0)
        tao1(ii)=LQSolvetao1(N,hn0(ii),0.025);
        tao2(ii)=LQSolvetao1(N,hn0(ii),0.05);
    end
    
    hp=figure;
    set(gcf,'Position',[714 301 616 486]);
    set(gcf,'PaperPositionMode','auto');
    h1=axes('position',[0.1150 0.0956 0.7750 0.8250]);
    axis equal
    xlim([0 1])
    ylim([0 1])
    hold on;
    box off;
    set(gca,'fontsize',FSa,'fontname',FNa);
    obs1=plot([0 1],[1 0],'k--','linewidth',LW1,'parent',h1);%gain=1�ο���
    obs2=plot(tao1,vn0,'r-.','linewidth',LW1,'parent',h1);%a=2.5%�ο���
    %plot(tao2,vn0,'c-.','linewidth',LW1,'parent',h1);
    obs3=plot([0,Coverrate],[1,Missrate],'b-','linewidth',LW3,'parent',h1);%����������
    obs4=scatter(Coverrate(indRR),Missrate(indRR),40,'o','markerfacecolor','m','markeredgecolor','k','linewidth',0.2);%Rֵ����
    obs5=scatter(Coverrate(indDists),Missrate(indDists),40,'o','markerfacecolor','g','markeredgecolor','k','linewidth',0.2);%����ԭ�������
    obs6=plot([0 1/Gains(indRR)],[1 0],'m--','linewidth',LW2,'parent',h1);
    obs7=plot([0 1/Gains(indDists)],[1 0],'g--','linewidth',LW2,'parent',h1);
    tmpstr1=['P1 Gain=',num2str(Gains(indRR),'%.4f')];
    tmpstr2=['P2 Gain=',num2str(Gains(indDists),'%.4f')];
    legend([obs1,obs2,obs3,obs4,obs5,obs6,obs7],{'Gain=1�ο���','\alpha=2.5%�ο���','����������','Rֵ����P1','����ԭ�������P2',tmpstr1,tmpstr2},'fontname',FNl);
    bb=Missrate(indRR)-Coverrate(indRR);
    %plot([Coverrate(indRR) (1-bb)/2],[Missrate(indRR) (1+bb)/2],'m:','linewidth',LW2,'parent',h1);
    %plot([0 Coverrate(indDists)],[0 Missrate(indDists)],'g:','linewidth',LW2,'parent',h1);
    text(0.4,0.5,['AUC=',num2str(AreaM)],'fontsize',FSa,'fontname',FNa);
    xlabel('������Ԫ��ʱ��ռ����\tau','fontsize',FSl,'fontname',FNl);
    ylabel('©����\nu','fontsize',FSl,'fontname',FNl);
    h2=axes('Position',get(gca,'Position'));
    set(h2,'color','none');
    axis equal
    xlim([0 1])
    ylim([0 1])
    set(gca,'fontsize',FSa,'fontname',FNa);
    set(h2,'yaxislocation','right');
    set(h2,'xaxislocation','top');
    set(h2,'xtick',[]);
    ntick=length(get(h2,'ytick'));
    dtick=ceil(length(h)/ntick);
    utick=h(end):-dtick:h(1);
    uutick=1-utick/N;
    set(h2,'ytick',uutick,'yticklabel',utick);
    ylabel('��׼��h','fontsize',FSl,'fontname',FNl);
    hold off;
    title({['����',num2str(ZhiBiao(mm,1)),'�����ϵ���, ���������ֵ:',num2str(ZhiBiao(mm,3)),', ���������ֵ:',num2str(ZhiBiao(mm,4)),', Ԥ��ʱ��:',num2str(ZhiBiao(mm,5)),'��, Ԥ����Χ:',num2str(ZhiBiao(mm,6)),'-',num2str(ZhiBiao(mm,7)),'km'];...
        ['Rֵ:',num2str(ZhiBiao(mm,8),'%.4f'),', R0ֵ:',num2str(ZhiBiao(mm,13),'%.4f'),', ʵ�ʷ���������:',num2str(ZhiBiao(mm,9)),', ��׼��:',num2str(ZhiBiao(mm,10),'%.4f'),', ��������:',num2str(ZhiBiao(mm,11),'%.4f'),', ����ˮƽ:',num2str(ZhiBiao(mm,12))]},...
        'FontName',FNl,'FontSize',FSl);
    Figname=[FnameQZ,'_MSI',num2str(ZhiBiao(mm,1))];%����ָ��ͼ��ÿ������ÿ���𼶵�һ��ͼ
    print(hp,[Figname,'.png'],'-dpng','-r600');
    saveas(hp,Figname,'pdf');
    close all;
end
end