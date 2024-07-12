%中国地震局地震预测研究所刘琦编制，最后调试时间2022-3-4，liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Molchan Error Diagram计算
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Threshold1,Threshold2,N,Hits,Coverrate,AreaM]=LQMolchanCompute2_ParNew(Utime,Udata,DTwindow,NoThreshold,UUseddzDateNum)
N=length(UUseddzDateNum);%总地震数
tmpThreshold=linspace(max(Udata),min(Udata),NoThreshold);%阈值网格线
Rlen=NoThreshold*(NoThreshold-1)/2;
ii1=NaN(1,Rlen);
ii2=NaN(1,Rlen);
kk=1;
for ik1=1:1:NoThreshold%负方向循环
    for ik2=1:1:ik1%正方向循环
        if ik1==ik2%始终位于Molchan图右下角点
            continue;
        end
        ii1(kk)=ik1;
        ii2(kk)=ik2;
        kk=kk+1;
    end
end

tmpThreshold1=NaN(1,length(ii1));
tmpThreshold2=NaN(1,length(ii1));
tmpCoverrate=NaN(1,length(ii1));
tmpHits=NaN(1,length(ii1));

parfor iUU=1:length(ii1)
    TT2=tmpThreshold(ii2(iUU));
    ind2=Udata>TT2;%大于正阈值
    TT1=tmpThreshold(ii1(iUU));
    ind1=Udata<TT1;%小于负阈值
    ind=ind1|ind2;%阈值组合
    CoUtime=Utime(ind);
    if isempty(CoUtime)
        continue;
    else
        AllCoUtime=repmat(CoUtime,1,DTwindow)+repmat([1:1:DTwindow],length(CoUtime),1);
        AllCoUtime=intersect(unique(AllCoUtime),Utime);
        tmp1=length(AllCoUtime)/(Utime(end)-Utime(1));%占有率
        tmp2=sum(ismember(UUseddzDateNum,AllCoUtime));%报准数
        
        tmpThreshold1(iUU)=TT1;
        tmpThreshold2(iUU)=TT2;
        tmpCoverrate(iUU)=tmp1;
        tmpHits(iUU)=tmp2;
    end
end

CM=[tmpCoverrate;tmpHits]';
[~,indUU]=unique(CM,'rows');%去重复点
tmpThreshold1=tmpThreshold1(indUU);
tmpThreshold2=tmpThreshold2(indUU);
tmpCoverrate=tmpCoverrate(indUU);
tmpHits=tmpHits(indUU);
[~,indUU]=sort(tmpCoverrate);%按时空占用率排序
tmpThreshold1=tmpThreshold1(indUU);
tmpThreshold2=tmpThreshold2(indUU);
tmpCoverrate=tmpCoverrate(indUU);
tmpHits=tmpHits(indUU);
[~,indUU]=unique(tmpHits);%挑各报准率对应的最优结果
Threshold1=tmpThreshold1(indUU);
Threshold2=tmpThreshold2(indUU);
Coverrate=tmpCoverrate(indUU);
Hits=tmpHits(indUU);
indUU=isnan(Hits)|Hits==0;
Threshold1(indUU)=[];%负方向阈值,即<Threshold1算异常
Threshold2(indUU)=[];%正方向阈值,即>Threshold2算异常
Coverrate(indUU)=[];%时间占有率，越小越好
Hits(indUU)=[];%报准地震数，越大越好

XX=[0,0,Coverrate,0];
YY=[0,1,1-Hits/N,0];
AreaM=polyarea(XX,YY);%曲线左下方面积，要小于0.5且越小越好
end