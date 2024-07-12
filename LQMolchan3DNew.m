%�й�����ֵ���Ԥ���о����������ƣ�������ʱ��2022-3-5��liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ü���Molchan���㺯�������ظ����ռ䴰���𼶴���ʱ�䴰���������������ŵ�ָ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BestThreshold1,BestThreshold2,BestRR,BestR0,BestEnum,BestHits,BestHitsrate,BestMiss,BestMissrate,BestCoverrate,BestGains,BestSig,BestDists,BestAreaM,SelectEC,NaNflag]=LQMolchan3DNew(DisRange,MM,DTwindow,NoThreshold,jwd,dzDate,dzLatitude,dzLongitude,dzMagnitude,dzDepth,dzLocation,Utime,Udata)
BestThreshold1=NaN(size(DisRange,1),size(MM,1),length(DTwindow));%����Ŀ¼ɸѡ����������µ�������ֵ����������ֵ
BestThreshold2=BestThreshold1;%����Ŀ¼ɸѡ����������µ�������ֵ����������ֵ
BestRR=NaN(size(DisRange,1),size(MM,1),length(DTwindow));%��ӦRֵ
BestR0=BestRR;%��ӦR0ֵ
BestEnum=BestRR;%��Ӧ��������
BestHits=BestRR;%��Ӧ��׼������
BestHitsrate=BestRR;%��Ӧ��׼��
BestMiss=BestRR;%��Ӧ©��������
BestMissrate=BestRR;%��Ӧ©����
BestCoverrate=BestRR;%��Ӧʱ��ռ����
BestGains=BestRR;%��Ӧ��������
BestSig=BestRR;%��Ӧ����ˮƽ
BestDists=BestRR;%��Ӧ����ԭ�����
BestAreaM=BestRR;%��Ӧ�������·�Χ�����
%��ά�±�ֱ��Ӧ�ռ䴰���𼶴���ʱ�䴰
%SelectECΪ���ݹ۲�ʱ����ɸѡ���ĵ���Ŀ¼����ά�±�ֱ��Ӧ�ռ䴰���𼶴�

NaNflag=0;%��ʾ��Ϊɸѡ����������޷�����
DZZS=0;%��ʾ����ʱ�䡢�ռ䡢��ɨ���ɸѡ���ĵ�������
%�۲�����ʱ�䷶Χ�ڵĵ���
[~,UseddzDate,~,UseddzLatitude,UseddzLongitude,UseddzMagnitude,UseddzDepth,UseddzLocation]=LQdzSelect1(jwd(2),jwd(1),dzDate,dzLatitude,dzLongitude,dzMagnitude,dzDepth,dzLocation,[],[],[Utime(1),Utime(end)]);
if isempty(UseddzMagnitude)%�۲�ʱ���޵���
    NaNflag=1;
    return;
end
for jj=1:1:size(DisRange,1)%�ռ�
    for kk=1:1:size(MM,1)%��
        %ɸѡ�ĵ���Ŀ¼
        [~,UUseddzDate,UUseddzDateNum,UUseddzLatitude,UUseddzLongitude,UUseddzMagnitude,UUseddzDepth,UUseddzLocation]=LQdzSelect1(jwd(2),jwd(1),UseddzDate,UseddzLatitude,UseddzLongitude,UseddzMagnitude,UseddzDepth,UseddzLocation,DisRange(jj,:),MM(kk,:),[]);
        SelectEC(jj,kk).Date=UUseddzDate;
        SelectEC(jj,kk).DateNum=UUseddzDateNum;
        SelectEC(jj,kk).Latitude=UUseddzLatitude;
        SelectEC(jj,kk).Longitude=UUseddzLongitude;
        SelectEC(jj,kk).Magnitude=UUseddzMagnitude;
        SelectEC(jj,kk).Depth=UUseddzDepth;
        SelectEC(jj,kk).Location=UUseddzLocation;
        SDZS=length(UUseddzMagnitude);
        DZZS=DZZS+SDZS;
        
        disp(['���ڽ��е�',num2str(jj),'-',num2str(size(DisRange,1)),'���ռ䴰����',num2str(kk),'-',num2str(size(MM,1)),'���𼶵�ɨ����㡣']);
        
        if SDZS==0%ɸѡ��������
            continue;
        end
        
        if kk>1&SDZS==length(SelectEC(jj,kk-1).Magnitude)%��һ���ռ䷶Χ�����е����Ǹ��𼶣�����������ɸѡ������Ŀ¼һ�£���Ӧָ��������䣬ʵ��Ӧ��Ӧ������Ԥ�����𼶵���
            BestThreshold1(jj,kk,:)=BestThreshold1(jj,kk-1,:);
            BestThreshold2(jj,kk,:)=BestThreshold2(jj,kk-1,:);
            BestHits(jj,kk,:)=BestHits(jj,kk-1,:);
            BestCoverrate(jj,kk,:)=BestCoverrate(jj,kk-1,:);
            BestSig(jj,kk,:)=BestSig(jj,kk-1,:);
            BestR0(jj,kk,:)=BestR0(jj,kk-1,:);
            continue;
        end
        
        if jj>1&SDZS==length(SelectEC(jj-1,kk).Magnitude)%�����Ȧ�㼯�У����Լ����������ɸѡ���ĵ��𲻱䣬���˿ռ�ռ���ʣ���Ӧָ�����������,ʵ��Ӧ����ѡ��С�ռ�
            BestThreshold1(jj,kk,:)=BestThreshold1(jj-1,kk,:);
            BestThreshold2(jj,kk,:)=BestThreshold2(jj-1,kk,:);
            BestHits(jj,kk,:)=BestHits(jj-1,kk,:);
            BestCoverrate(jj,kk,:)=BestCoverrate(jj-1,kk,:);
            BestSig(jj,kk,:)=BestSig(jj-1,kk,:);
            BestR0(jj,kk,:)=BestR0(jj-1,kk,:);
            continue;      
        end
        
        for ll=1:1:length(DTwindow)%ʱ�䴰
            disp([num2str(ll),'-',num2str(length(DTwindow)),'ʱ�䴰']);
            [Threshold1,Threshold2,N,Hits,Coverrate,AreaM]=LQMolchanCompute2_ParNew(Utime,Udata,DTwindow(ll),NoThreshold,UUseddzDateNum);
            BestEnum(jj,kk,ll)=N;
            BestAreaM(jj,kk,ll)=AreaM;
            RR=Hits/N-Coverrate;
            Dists=sqrt((1-Hits/N).^2+Coverrate.^2);
            [~,indRR]=max(RR);
            [~,indDists]=min(Dists);
            %�����ռ䡢�𼶡�ʱ��ȴ���������µ�������ֵѡȡ��
            indAU=union(indRR,indDists);%��ѡRֵ�������ԭ����������е�
            if length(indAU)==1%ֻ��1����
                BestThreshold1(jj,kk,ll)=Threshold1(indAU);
                BestThreshold2(jj,kk,ll)=Threshold2(indAU);
                BestHits(jj,kk,ll)=Hits(indAU);
                BestCoverrate(jj,kk,ll)=Coverrate(indAU);
                BestSig(jj,kk,ll)=LQSignificanceLevelT(N,Hits(indAU),Coverrate(indAU),2);
                BestR0(jj,kk,ll)=LQGetR00(Hits(indAU),N-Hits(indAU));
            elseif length(indAU)>1%�����Ψһ�����ո������������ˮƽ�͡�R0С��ԭ������ɸѡ��
                %������ͬ�𼶷�Χ���ռ䷶Χ�������������ͬһ�׵�������������£�������Rֵ��ͬ�Ҹ���������ͬ����������ˮƽ��R0ֵ�ض���ͬ��
                tmpThreshold1=Threshold1(indAU);
                tmpThreshold2=Threshold2(indAU);
                tmpHits=Hits(indAU);
                tmpCoverrate=Coverrate(indAU);
                tmpRR=RR(indAU);
                tmpGain=tmpHits/N./tmpCoverrate;
                for nnn=1:1:length(indAU)
                    tmpSig(nnn)=LQSignificanceLevelT(N,tmpHits(nnn),tmpCoverrate(nnn),2);
                    tmpR0(nnn)=LQGetR00(tmpHits(nnn),N-tmpHits(nnn));
                end
                [~,indAU]=max(tmpGain);%�����������
                if length(indAU)==1%Ψһ
                    BestThreshold1(jj,kk,ll)=tmpThreshold1(indAU);
                    BestThreshold2(jj,kk,ll)=tmpThreshold2(indAU);
                    BestHits(jj,kk,ll)=tmpHits(indAU);
                    BestCoverrate(jj,kk,ll)=tmpCoverrate(indAU);
                    BestSig(jj,kk,ll)=tmpSig(indAU);
                    BestR0(jj,kk,ll)=tmpR0(indAU);
                else%��Ψһ����ʵ���������ٱȽ���
                    tmpThreshold1=tmpThreshold1(indAU);
                    tmpThreshold2=tmpThreshold2(indAU);
                    tmpHits=tmpHits(indAU);
                    tmpCoverrate=tmpCoverrate(indAU);
                    tmpSig=tmpSig(indAU);
                    tmpR0=tmpR0(indAU);
                    [~,indAU]=min(tmpSig);%������ˮƽ���
                    if length(indAU)==1%Ψһ
                        BestThreshold1(jj,kk,ll)=tmpThreshold1(indAU);
                        BestThreshold2(jj,kk,ll)=tmpThreshold2(indAU);
                        BestHits(jj,kk,ll)=tmpHits(indAU);
                        BestCoverrate(jj,kk,ll)=tmpCoverrate(indAU);
                        BestSig(jj,kk,ll)=tmpSig(indAU);
                        BestR0(jj,kk,ll)=tmpR0(indAU);
                    else%��Ψһ
                        tmpThreshold1=tmpThreshold1(indAU);
                        tmpThreshold2=tmpThreshold2(indAU);
                        tmpHits=tmpHits(indAU);
                        tmpCoverrate=tmpCoverrate(indAU);
                        tmpSig=tmpSig(indAU);
                        tmpR0=tmpR0(indAU);
                        [~,indAU]=min(tmpR0);%R0��С
                        BestThreshold1(jj,kk,ll)=tmpThreshold1(indAU(1));
                        BestThreshold2(jj,kk,ll)=tmpThreshold2(indAU(1));
                        BestHits(jj,kk,ll)=tmpHits(indAU(1));
                        BestCoverrate(jj,kk,ll)=tmpCoverrate(indAU(1));
                        BestSig(jj,kk,ll)=tmpSig(indAU(1));
                        BestR0(jj,kk,ll)=tmpR0(indAU(1));
                    end
                end
            else
                continue;
            end
        end
    end
end
BestHitsrate=BestHits./BestEnum;
BestRR=BestHitsrate-BestCoverrate;
BestMiss=BestEnum-BestHits;
BestMissrate=1-BestHitsrate;
BestGains=BestHitsrate./BestCoverrate;
BestDists=sqrt(BestMissrate.^2+BestCoverrate.^2);

if DZZS==0%�޵���
    NaNflag=1;
end
end