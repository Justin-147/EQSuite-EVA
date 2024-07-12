%�й�����ֵ���Ԥ���о����������ƣ�������ʱ��2022-3-4��liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Molchan Error Diagram����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Threshold1,Threshold2,N,Hits,Coverrate,AreaM]=LQMolchanCompute2_ParNew(Utime,Udata,DTwindow,NoThreshold,UUseddzDateNum)
N=length(UUseddzDateNum);%�ܵ�����
tmpThreshold=linspace(max(Udata),min(Udata),NoThreshold);%��ֵ������
Rlen=NoThreshold*(NoThreshold-1)/2;
ii1=NaN(1,Rlen);
ii2=NaN(1,Rlen);
kk=1;
for ik1=1:1:NoThreshold%������ѭ��
    for ik2=1:1:ik1%������ѭ��
        if ik1==ik2%ʼ��λ��Molchanͼ���½ǵ�
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
    ind2=Udata>TT2;%��������ֵ
    TT1=tmpThreshold(ii1(iUU));
    ind1=Udata<TT1;%С�ڸ���ֵ
    ind=ind1|ind2;%��ֵ���
    CoUtime=Utime(ind);
    if isempty(CoUtime)
        continue;
    else
        AllCoUtime=repmat(CoUtime,1,DTwindow)+repmat([1:1:DTwindow],length(CoUtime),1);
        AllCoUtime=intersect(unique(AllCoUtime),Utime);
        tmp1=length(AllCoUtime)/(Utime(end)-Utime(1));%ռ����
        tmp2=sum(ismember(UUseddzDateNum,AllCoUtime));%��׼��
        
        tmpThreshold1(iUU)=TT1;
        tmpThreshold2(iUU)=TT2;
        tmpCoverrate(iUU)=tmp1;
        tmpHits(iUU)=tmp2;
    end
end

CM=[tmpCoverrate;tmpHits]';
[~,indUU]=unique(CM,'rows');%ȥ�ظ���
tmpThreshold1=tmpThreshold1(indUU);
tmpThreshold2=tmpThreshold2(indUU);
tmpCoverrate=tmpCoverrate(indUU);
tmpHits=tmpHits(indUU);
[~,indUU]=sort(tmpCoverrate);%��ʱ��ռ��������
tmpThreshold1=tmpThreshold1(indUU);
tmpThreshold2=tmpThreshold2(indUU);
tmpCoverrate=tmpCoverrate(indUU);
tmpHits=tmpHits(indUU);
[~,indUU]=unique(tmpHits);%������׼�ʶ�Ӧ�����Ž��
Threshold1=tmpThreshold1(indUU);
Threshold2=tmpThreshold2(indUU);
Coverrate=tmpCoverrate(indUU);
Hits=tmpHits(indUU);
indUU=isnan(Hits)|Hits==0;
Threshold1(indUU)=[];%��������ֵ,��<Threshold1���쳣
Threshold2(indUU)=[];%��������ֵ,��>Threshold2���쳣
Coverrate(indUU)=[];%ʱ��ռ���ʣ�ԽСԽ��
Hits(indUU)=[];%��׼��������Խ��Խ��

XX=[0,0,Coverrate,0];
YY=[0,1,1-Hits/N,0];
AreaM=polyarea(XX,YY);%�������·������ҪС��0.5��ԽСԽ��
end