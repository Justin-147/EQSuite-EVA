%中国地震局地震预测研究所刘琦编制，最后调试时间2022-3-5，liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%调用计算Molchan计算函数，返回给定空间窗、震级窗、时间窗参数组合情况下最优的指标量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BestThreshold1,BestThreshold2,BestRR,BestR0,BestEnum,BestHits,BestHitsrate,BestMiss,BestMissrate,BestCoverrate,BestGains,BestSig,BestDists,BestAreaM,SelectEC,NaNflag]=LQMolchan3DNew(DisRange,MM,DTwindow,NoThreshold,jwd,dzDate,dzLatitude,dzLongitude,dzMagnitude,dzDepth,dzLocation,Utime,Udata)
BestThreshold1=NaN(size(DisRange,1),size(MM,1),length(DTwindow));%地震目录筛选参数各组合下的最优阈值，负方向阈值
BestThreshold2=BestThreshold1;%地震目录筛选参数各组合下的最优阈值，正方向阈值
BestRR=NaN(size(DisRange,1),size(MM,1),length(DTwindow));%对应R值
BestR0=BestRR;%对应R0值
BestEnum=BestRR;%对应地震总数
BestHits=BestRR;%对应报准地震数
BestHitsrate=BestRR;%对应报准率
BestMiss=BestRR;%对应漏报地震数
BestMissrate=BestRR;%对应漏报率
BestCoverrate=BestRR;%对应时空占用率
BestGains=BestRR;%对应概率增益
BestSig=BestRR;%对应显著水平
BestDists=BestRR;%对应距离原点距离
BestAreaM=BestRR;%对应曲线左下方围限面积
%三维下标分别对应空间窗、震级窗、时间窗
%SelectEC为数据观测时段内筛选出的地震目录，二维下标分别对应空间窗和震级窗

NaNflag=0;%表示因为筛选不到地震而无法评估
DZZS=0;%表示所有时间、空间、震级扫描后，筛选到的地震总数
%观测数据时间范围内的地震
[~,UseddzDate,~,UseddzLatitude,UseddzLongitude,UseddzMagnitude,UseddzDepth,UseddzLocation]=LQdzSelect1(jwd(2),jwd(1),dzDate,dzLatitude,dzLongitude,dzMagnitude,dzDepth,dzLocation,[],[],[Utime(1),Utime(end)]);
if isempty(UseddzMagnitude)%观测时段无地震
    NaNflag=1;
    return;
end
for jj=1:1:size(DisRange,1)%空间
    for kk=1:1:size(MM,1)%震级
        %筛选的地震目录
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
        
        disp(['正在进行第',num2str(jj),'-',num2str(size(DisRange,1)),'个空间窗，第',num2str(kk),'-',num2str(size(MM,1)),'个震级的扫描计算。']);
        
        if SDZS==0%筛选不出地震
            continue;
        end
        
        if kk>1&SDZS==length(SelectEC(jj,kk-1).Magnitude)%在一定空间范围内所有地震都是高震级，所以两档震级筛选出来的目录一致，对应指标参数不变，实际应用应该用来预报高震级地震
            BestThreshold1(jj,kk,:)=BestThreshold1(jj,kk-1,:);
            BestThreshold2(jj,kk,:)=BestThreshold2(jj,kk-1,:);
            BestHits(jj,kk,:)=BestHits(jj,kk-1,:);
            BestCoverrate(jj,kk,:)=BestCoverrate(jj,kk-1,:);
            BestSig(jj,kk,:)=BestSig(jj,kk-1,:);
            BestR0(jj,kk,:)=BestR0(jj,kk-1,:);
            continue;
        end
        
        if jj>1&SDZS==length(SelectEC(jj-1,kk).Magnitude)%地震分圈层集中，所以继续扩大面积筛选到的地震不变，除了空间占有率，对应指标参数均不变,实际应用中选最小空间
            BestThreshold1(jj,kk,:)=BestThreshold1(jj-1,kk,:);
            BestThreshold2(jj,kk,:)=BestThreshold2(jj-1,kk,:);
            BestHits(jj,kk,:)=BestHits(jj-1,kk,:);
            BestCoverrate(jj,kk,:)=BestCoverrate(jj-1,kk,:);
            BestSig(jj,kk,:)=BestSig(jj-1,kk,:);
            BestR0(jj,kk,:)=BestR0(jj-1,kk,:);
            continue;      
        end
        
        for ll=1:1:length(DTwindow)%时间窗
            disp([num2str(ll),'-',num2str(length(DTwindow)),'时间窗']);
            [Threshold1,Threshold2,N,Hits,Coverrate,AreaM]=LQMolchanCompute2_ParNew(Utime,Udata,DTwindow(ll),NoThreshold,UUseddzDateNum);
            BestEnum(jj,kk,ll)=N;
            BestAreaM(jj,kk,ll)=AreaM;
            RR=Hits/N-Coverrate;
            Dists=sqrt((1-Hits/N).^2+Coverrate.^2);
            [~,indRR]=max(RR);
            [~,indDists]=min(Dists);
            %给定空间、震级、时间等窗参数组合下的最优阈值选取。
            indAU=union(indRR,indDists);%挑选R值最大或距离原点最近的所有点
            if length(indAU)==1%只有1个点
                BestThreshold1(jj,kk,ll)=Threshold1(indAU);
                BestThreshold2(jj,kk,ll)=Threshold2(indAU);
                BestHits(jj,kk,ll)=Hits(indAU);
                BestCoverrate(jj,kk,ll)=Coverrate(indAU);
                BestSig(jj,kk,ll)=LQSignificanceLevelT(N,Hits(indAU),Coverrate(indAU),2);
                BestR0(jj,kk,ll)=LQGetR00(Hits(indAU),N-Hits(indAU));
            elseif length(indAU)>1%如果不唯一，则按照概率增益大、显著水平低、R0小的原则依次筛选。
                %对于相同震级范围、空间范围的情况，即采用同一套地震样本的情况下，若两点R值相同且概率增益相同，则其显著水平、R0值必定相同。
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
                [~,indAU]=max(tmpGain);%概率增益最大
                if length(indAU)==1%唯一
                    BestThreshold1(jj,kk,ll)=tmpThreshold1(indAU);
                    BestThreshold2(jj,kk,ll)=tmpThreshold2(indAU);
                    BestHits(jj,kk,ll)=tmpHits(indAU);
                    BestCoverrate(jj,kk,ll)=tmpCoverrate(indAU);
                    BestSig(jj,kk,ll)=tmpSig(indAU);
                    BestR0(jj,kk,ll)=tmpR0(indAU);
                else%不唯一，其实后续不用再比较了
                    tmpThreshold1=tmpThreshold1(indAU);
                    tmpThreshold2=tmpThreshold2(indAU);
                    tmpHits=tmpHits(indAU);
                    tmpCoverrate=tmpCoverrate(indAU);
                    tmpSig=tmpSig(indAU);
                    tmpR0=tmpR0(indAU);
                    [~,indAU]=min(tmpSig);%显著性水平最低
                    if length(indAU)==1%唯一
                        BestThreshold1(jj,kk,ll)=tmpThreshold1(indAU);
                        BestThreshold2(jj,kk,ll)=tmpThreshold2(indAU);
                        BestHits(jj,kk,ll)=tmpHits(indAU);
                        BestCoverrate(jj,kk,ll)=tmpCoverrate(indAU);
                        BestSig(jj,kk,ll)=tmpSig(indAU);
                        BestR0(jj,kk,ll)=tmpR0(indAU);
                    else%不唯一
                        tmpThreshold1=tmpThreshold1(indAU);
                        tmpThreshold2=tmpThreshold2(indAU);
                        tmpHits=tmpHits(indAU);
                        tmpCoverrate=tmpCoverrate(indAU);
                        tmpSig=tmpSig(indAU);
                        tmpR0=tmpR0(indAU);
                        [~,indAU]=min(tmpR0);%R0最小
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

if DZZS==0%无地震
    NaNflag=1;
end
end