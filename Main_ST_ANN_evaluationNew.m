%基于S变换提取破年变信息，并结合R值和Molchan图表对预测指标进行效能评估
%作者：刘琦
%单位：中国地震局地震预测研究所
%准备好原始数据、地震目录，修改索引文件，然后运行
%该版本未考虑跨断层数据、干扰信息表、给定多边形地震筛选边界、自动绘制空间图
%Extracting the annual cycle breaking information based on S
%transformation,and evaluating the performance of the predicition indicator
%by the combination of the R value and the Molchan diagram
%Author: Liu Qi
%Institution: IEF
%Contact: liuqi@ief.ac.cn
%Last Modified: 2022/3/5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%基本参数
clear;clc;close all;
EarthquakeFile='全国地震目录Ms5_wp.EQT';%所用地震目录文件名
FigPathU='Fig';
ResultPathU='Result';
mkdir(FigPathU);%图件存放路径
mkdir(ResultPathU);%结果存放路径
%年变提取参数
QS=999999;%缺数标记
ann_cyc2=[500,14];%计算周期带
ann_cyc1=[500,250];%年变周期带
minLen=2;%数据最小长度，年
SPadLen=0.075;%S变换前单侧扩边数比例
%地震目录筛选参数，按照距离、震级、时窗等筛选目录
DisRange=50:50:500;%筛选的最大距离，km
DisRange=[zeros(length(DisRange),1),DisRange'];%筛选的距离范围，目前为实心圆范围，也可更改为圆环范围
MM=5:1:7;%筛选的最小震级
MM=[MM',Inf(length(MM),1)];%筛选的震级范围
DTwindow=[7,14,30:30:1080];%时窗范围，单位天
NoThreshold=100;%阈值扫描个数,至少3个
%绘图输出等参数
InMCIMode=1;%全参数空间扫描评估模式，1为计算，其它设置为读取已有结果文件
PlotLanMode=1;%图件标注的格式，1为中文，其它设置为英文
FlagPanC=2;%年变显著度检查图，1全频段绘制，2截断频段绘制，其它设置不绘制
FlagPanE=1;%原始观测及ANA、ONA曲线图，1绘制，其它设置不绘制
FlagOanE=1;%ANA及ONA时序文本数据输出，1输出，其它设置不输出
FlagPanAP=1;%综合指标图，1绘制，其它设置不绘制
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Olon,Olat,OFile]=textread('索引文件.txt','%f %f %s');%测项经度、纬度、文件名
NFZ=length(OFile);%NFZ为待处理文件个数
if NFZ==0%索引文件中无内容
    return;
end

number_null=0;%无效数据文件数量
number_downs=0;%采样率不够的文件数量
number_slen=0;%数据长度过短的文件数量
number_upc=0;%扩边或截边未成功文件数量
number_uann=0;%年变不清晰文件数量

[eqDateTime,eqLat,eqLon,eqMag,eqDepth,eqLocation]=ReadEQT_lq(EarthquakeFile);%s,f,f,f,f,s

for iiNFZ=1:1:NFZ
    dbfile=OFile{iiNFZ};
    [Ptmp,Ftmp,Etmp]=fileparts(dbfile);
    FF=[Ftmp,Etmp];
    disp(['文件总数:',num2str(NFZ),'，目前正在处理第',num2str(iiNFZ),'个，文件名为：',FF]);
    
    %导入数据并转换成标准2列格式
    [time_s,data_s]=Data_Trans_lq(dbfile);
    
    %数据预处理：补断数、去台阶和突跳、邻近值简单进行缺数插值
    [time_p,data_p,oflag,wz_QS]=Data_Prep_lq(time_s,data_s,QS);
    if oflag==0%全部数据无效
        number_null=number_null+1;
        file_null(number_null)={dbfile};
        continue;
    end
    time_QS=time_p(wz_QS);
    
    %数据降采样，降为日值
    [time_p,data_p,oflag]=Data_DownS_lq(time_p,data_p,'整日值');
    if oflag==0%全部数据无效
        number_downs=number_downs+1;
        file_downs(number_downs)={dbfile};
        continue;
    end
    len1=length(num2str(time_p(1)));
    len2=length(num2str(time_QS(1)));
    time_QS_D=time_QS(mod(time_QS,10^(len2-len1))==0)/10^(len2-len1);
    [~,wz_QS,~]=intersect(time_p,time_QS_D);%日值里缺数位置
    
    %数据长度过短跳过
    if length(data_p)<minLen*365%数据长度过短
        number_slen=number_slen+1;
        file_slen(number_slen)={dbfile};
        continue;
    end
    
    %小波滤波,不需要
    %data_p=Data_WaveFilt(data_p);
    
    %年变不清晰跳过
    Ratio_Ann=Data_Ann_CheckP(data_p,ann_cyc1,FlagPanC,PlotLanMode,FF,[pwd,'\Fig\']);
    if Ratio_Ann<1
        number_uann=number_uann+1;
        file_uann(number_uann)={dbfile};
        continue;
    end
    
    %考虑边界效应，数据需要扩边
    width_pad=ceil(length(data_p)*SPadLen);%单侧扩边样点数
    type_pad='all';%双侧扩边，可选参数'all','left','right'
    flag_pad='pad';%'pad'扩边,'cut'截边
    [data_p,oflag]=Data_PadorCut(data_p,width_pad,type_pad,flag_pad);
    if oflag==0%扩边未成功
        number_upc=number_upc+1;
        file_upc(number_upc)={dbfile};
        continue;
    end
    
    %S变换
    [ST,~,Fre]=stm1(data_p,round(length(data_p)/ann_cyc2(1)),round(length(data_p)/ann_cyc2(2)));%ST结果为复数,周期范围为500-14天左右
    Cyc=1./Fre;%信号周期，单位同采样率单位
    [~,ind_ann_cyc1]=min(abs(repmat(Cyc',1,2)-repmat(ann_cyc1,length(Cyc),1)));%年变周期带索引，周期范围为500-250天左右
    
    %考虑边界效应，数据需要截边
    flag_pad='cut';%'pad'扩边,'cut'截边
    [data_p,~]=Data_PadorCut(data_p,width_pad,type_pad,flag_pad);
    [ST,~]=Data_PadorCut(ST',width_pad,type_pad,flag_pad);
    ST=ST';
    
    %整理结果
    [~,data_mean_ST_ann,data_mean_ST_other,~,~,~]=Data_out(data_p,ST,Cyc,Fre,ind_ann_cyc1);
    
    f_nn=find(FF=='.',1,'last')-1;
    if FlagPanE==1
        %绘图并保存
        data_p(wz_QS)=NaN;
        data_mean_ST_ann(wz_QS)=NaN;
        data_mean_ST_other(wz_QS)=NaN;
        Plot_Ann(time_p,data_p,data_mean_ST_ann,data_mean_ST_other,ann_cyc1,ann_cyc2,PlotLanMode);
        Figname=[[pwd,'\',FigPathU,'\'],FF(1:f_nn),'_破年变信息提取'];
        print(gcf,[Figname,'.png'],'-dpng','-r600');
        saveas(gcf,Figname,'pdf');
        close all;
    end
    
    if FlagOanE==1
        %输出
        data_mean_ST_ann(wz_QS)=QS;
        data_mean_ST_other(wz_QS)=QS;
        fm=['%8i %.5f\n'];
        outname=strcat(pwd,'\',ResultPathU,'\',FF(1:f_nn),'_ST_ANA','.txt');
        fidof=fopen(outname,'wt');
        fprintf(fidof,fm,[time_p';data_mean_ST_ann]);
        fclose(fidof);
        outname=strcat(pwd,'\',ResultPathU,'\',FF(1:f_nn),'_ST_ONA','.txt');
        fidof=fopen(outname,'wt');
        fprintf(fidof,fm,[time_p';data_mean_ST_other]);
        fclose(fidof);
    end
    
    StaLatLon=[Olon(iiNFZ),Olat(iiNFZ)];%台站经纬度
    Utime=datenum(num2str(time_p),'yyyymmdd');
    outname=strcat(pwd,'\',ResultPathU,'\',FF(1:f_nn),'_Basic_result.mat');
    save(outname,'ann_cyc1','ann_cyc2','Cyc','data_mean_ST_ann','data_mean_ST_other','data_p','DisRange','DTwindow','FF','Fre','ind_ann_cyc1','minLen','MM','NoThreshold','QS','Ratio_Ann','SPadLen','ST','StaLatLon','time_p','wz_QS');
    
    outApp='Eval_other';
    outRt=strcat(pwd,'\',ResultPathU,'\',FF(1:f_nn),'_',outApp);
    outFt=strcat(pwd,'\',FigPathU,'\',FF(1:f_nn),'_',outApp);
    if InMCIMode==1
        Udata=data_mean_ST_other;
        Udata(wz_QS)=NaN;
        tt1=tic;
        [BestThreshold1,BestThreshold2,BestRR,BestR0,BestEnum,BestHits,BestHitsrate,BestMiss,BestMissrate,BestCoverrate,BestGains,BestSig,BestDists,BestAreaM,SelectEC,NaNflag]=LQMolchan3DNew(DisRange,MM,DTwindow,NoThreshold,StaLatLon,str2num(eqDateTime),eqLat,eqLon,eqMag,eqDepth,eqLocation,Utime,Udata);
        ttelapsed1=toc(tt1);
        if NaNflag==1%无地震
            disp([outApp,': No earthquakes!']);
        else
            disp([outApp,' uses ',num2str(ttelapsed1),'seconds.']);
            outname=[outRt,'.mat'];
            save(outname,'BestThreshold1','BestThreshold2','BestRR','BestR0','BestEnum','BestHits','BestHitsrate','BestMiss','BestMissrate','BestCoverrate','BestGains','BestSig','BestDists','BestAreaM','SelectEC','NaNflag');
        end
    else
        outname=[outRt,'.mat'];
        load(outname);
        Udata=data_mean_ST_other;
        Udata(wz_QS)=NaN;
    end
    
    outname=strcat(pwd,'\',ResultPathU,'\',FF(1:f_nn),'_BestIndex.mat');
    load(outname);
    ZhiBiao1=LQZhiBiaoComp1(BestRR,BestThreshold1,BestThreshold2,BestHitsrate,BestR0,BestGains,BestSig,BestEnum,DisRange,DTwindow,MM,outFt,FlagPanAP);    
    YL={'其它频段归一化平均幅度ONA';['(周期',num2str(ann_cyc2(2)),'-',num2str(ann_cyc2(1)),'天)']};
    OutDZML1=LQPlotSingleCurve1(ZhiBiao1,YL,SelectEC,DisRange,Utime,Udata,outFt);
    LQMolchan1DNew(ZhiBiao1,SelectEC,Utime,Udata,NoThreshold,DisRange,outFt);
    clear SelectEC
    
    outApp='Eval_ann';
    outRt=strcat(pwd,'\',ResultPathU,'\',FF(1:f_nn),'_',outApp);
    outFt=strcat(pwd,'\',FigPathU,'\',FF(1:f_nn),'_',outApp);
    if InMCIMode==1
        Udata=data_mean_ST_ann;
        Udata(wz_QS)=NaN;
        tt2=tic;
        [BestThreshold1,BestThreshold2,BestRR,BestR0,BestEnum,BestHits,BestHitsrate,BestMiss,BestMissrate,BestCoverrate,BestGains,BestSig,BestDists,BestAreaM,SelectEC,NaNflag]=LQMolchan3DNew(DisRange,MM,DTwindow,NoThreshold,StaLatLon,str2num(eqDateTime),eqLat,eqLon,eqMag,eqDepth,eqLocation,Utime,Udata);
        ttelapsed2=toc(tt2);
        if NaNflag==1%无地震
            disp([outApp,': No earthquakes!']);
        else
            disp([outApp,' uses ',num2str(ttelapsed2),'seconds.']);
            outname=[outRt,'.mat'];
            save(outname,'BestThreshold1','BestThreshold2','BestRR','BestR0','BestEnum','BestHits','BestHitsrate','BestMiss','BestMissrate','BestCoverrate','BestGains','BestSig','BestDists','BestAreaM','SelectEC','NaNflag');
        end
    else
        outname=[outRt,'.mat'];
        load(outname);
        Udata=data_mean_ST_ann;
        Udata(wz_QS)=NaN;
    end
    ZhiBiao2=LQZhiBiaoComp1(BestRR,BestThreshold1,BestThreshold2,BestHitsrate,BestR0,BestGains,BestSig,BestEnum,DisRange,DTwindow,MM,outFt,FlagPanAP);
    YL={'年频段归一化平均幅度ANA';['(周期',num2str(ann_cyc1(2)),'-',num2str(ann_cyc1(1)),'天)']};
    OutDZML2=LQPlotSingleCurve1(ZhiBiao2,YL,SelectEC,DisRange,Utime,Udata,outFt);
    LQMolchan1DNew(ZhiBiao2,SelectEC,Utime,Udata,NoThreshold,DisRange,outFt);
    clear SelectEC
    
    outname=strcat(pwd,'\',ResultPathU,'\',FF(1:f_nn),'_BestIndex.mat');
    save(outname,'ZhiBiao1','ZhiBiao2','OutDZML1','OutDZML2');
    clear OutDZML1 OutDZML2
end
%LQMolchanCompute2_ParNew采用了并行命令，最好结束前终止Matlab并行计算环境
delete(gcp('nocreate'));