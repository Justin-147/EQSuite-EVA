%����S�任��ȡ�������Ϣ�������Rֵ��Molchanͼ���Ԥ��ָ�����Ч������
%���ߣ�����
%��λ���й�����ֵ���Ԥ���о���
%׼����ԭʼ���ݡ�����Ŀ¼���޸������ļ���Ȼ������
%�ð汾δ���ǿ�ϲ����ݡ�������Ϣ����������ε���ɸѡ�߽硢�Զ����ƿռ�ͼ
%Extracting the annual cycle breaking information based on S
%transformation,and evaluating the performance of the predicition indicator
%by the combination of the R value and the Molchan diagram
%Author: Liu Qi
%Institution: IEF
%Contact: liuqi@ief.ac.cn
%Last Modified: 2022/3/5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��������
clear;clc;close all;
EarthquakeFile='ȫ������Ŀ¼Ms5_wp.EQT';%���õ���Ŀ¼�ļ���
FigPathU='Fig';
ResultPathU='Result';
mkdir(FigPathU);%ͼ�����·��
mkdir(ResultPathU);%������·��
%�����ȡ����
QS=999999;%ȱ�����
ann_cyc2=[500,14];%�������ڴ�
ann_cyc1=[500,250];%������ڴ�
minLen=2;%������С���ȣ���
SPadLen=0.075;%S�任ǰ��������������
%����Ŀ¼ɸѡ���������վ��롢�𼶡�ʱ����ɸѡĿ¼
DisRange=50:50:500;%ɸѡ�������룬km
DisRange=[zeros(length(DisRange),1),DisRange'];%ɸѡ�ľ��뷶Χ��ĿǰΪʵ��Բ��Χ��Ҳ�ɸ���ΪԲ����Χ
MM=5:1:7;%ɸѡ����С��
MM=[MM',Inf(length(MM),1)];%ɸѡ���𼶷�Χ
DTwindow=[7,14,30:30:1080];%ʱ����Χ����λ��
NoThreshold=100;%��ֵɨ�����,����3��
%��ͼ����Ȳ���
InMCIMode=1;%ȫ�����ռ�ɨ������ģʽ��1Ϊ���㣬��������Ϊ��ȡ���н���ļ�
PlotLanMode=1;%ͼ����ע�ĸ�ʽ��1Ϊ���ģ���������ΪӢ��
FlagPanC=2;%��������ȼ��ͼ��1ȫƵ�λ��ƣ�2�ض�Ƶ�λ��ƣ��������ò�����
FlagPanE=1;%ԭʼ�۲⼰ANA��ONA����ͼ��1���ƣ��������ò�����
FlagOanE=1;%ANA��ONAʱ���ı����������1������������ò����
FlagPanAP=1;%�ۺ�ָ��ͼ��1���ƣ��������ò�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Olon,Olat,OFile]=textread('�����ļ�.txt','%f %f %s');%����ȡ�γ�ȡ��ļ���
NFZ=length(OFile);%NFZΪ�������ļ�����
if NFZ==0%�����ļ���������
    return;
end

number_null=0;%��Ч�����ļ�����
number_downs=0;%�����ʲ������ļ�����
number_slen=0;%���ݳ��ȹ��̵��ļ�����
number_upc=0;%���߻�ر�δ�ɹ��ļ�����
number_uann=0;%��䲻�����ļ�����

[eqDateTime,eqLat,eqLon,eqMag,eqDepth,eqLocation]=ReadEQT_lq(EarthquakeFile);%s,f,f,f,f,s

for iiNFZ=1:1:NFZ
    dbfile=OFile{iiNFZ};
    [Ptmp,Ftmp,Etmp]=fileparts(dbfile);
    FF=[Ftmp,Etmp];
    disp(['�ļ�����:',num2str(NFZ),'��Ŀǰ���ڴ����',num2str(iiNFZ),'�����ļ���Ϊ��',FF]);
    
    %�������ݲ�ת���ɱ�׼2�и�ʽ
    [time_s,data_s]=Data_Trans_lq(dbfile);
    
    %����Ԥ������������ȥ̨�׺�ͻ�����ڽ�ֵ�򵥽���ȱ����ֵ
    [time_p,data_p,oflag,wz_QS]=Data_Prep_lq(time_s,data_s,QS);
    if oflag==0%ȫ��������Ч
        number_null=number_null+1;
        file_null(number_null)={dbfile};
        continue;
    end
    time_QS=time_p(wz_QS);
    
    %���ݽ���������Ϊ��ֵ
    [time_p,data_p,oflag]=Data_DownS_lq(time_p,data_p,'����ֵ');
    if oflag==0%ȫ��������Ч
        number_downs=number_downs+1;
        file_downs(number_downs)={dbfile};
        continue;
    end
    len1=length(num2str(time_p(1)));
    len2=length(num2str(time_QS(1)));
    time_QS_D=time_QS(mod(time_QS,10^(len2-len1))==0)/10^(len2-len1);
    [~,wz_QS,~]=intersect(time_p,time_QS_D);%��ֵ��ȱ��λ��
    
    %���ݳ��ȹ�������
    if length(data_p)<minLen*365%���ݳ��ȹ���
        number_slen=number_slen+1;
        file_slen(number_slen)={dbfile};
        continue;
    end
    
    %С���˲�,����Ҫ
    %data_p=Data_WaveFilt(data_p);
    
    %��䲻��������
    Ratio_Ann=Data_Ann_CheckP(data_p,ann_cyc1,FlagPanC,PlotLanMode,FF,[pwd,'\Fig\']);
    if Ratio_Ann<1
        number_uann=number_uann+1;
        file_uann(number_uann)={dbfile};
        continue;
    end
    
    %���Ǳ߽�ЧӦ��������Ҫ����
    width_pad=ceil(length(data_p)*SPadLen);%��������������
    type_pad='all';%˫�����ߣ���ѡ����'all','left','right'
    flag_pad='pad';%'pad'����,'cut'�ر�
    [data_p,oflag]=Data_PadorCut(data_p,width_pad,type_pad,flag_pad);
    if oflag==0%����δ�ɹ�
        number_upc=number_upc+1;
        file_upc(number_upc)={dbfile};
        continue;
    end
    
    %S�任
    [ST,~,Fre]=stm1(data_p,round(length(data_p)/ann_cyc2(1)),round(length(data_p)/ann_cyc2(2)));%ST���Ϊ����,���ڷ�ΧΪ500-14������
    Cyc=1./Fre;%�ź����ڣ���λͬ�����ʵ�λ
    [~,ind_ann_cyc1]=min(abs(repmat(Cyc',1,2)-repmat(ann_cyc1,length(Cyc),1)));%������ڴ����������ڷ�ΧΪ500-250������
    
    %���Ǳ߽�ЧӦ��������Ҫ�ر�
    flag_pad='cut';%'pad'����,'cut'�ر�
    [data_p,~]=Data_PadorCut(data_p,width_pad,type_pad,flag_pad);
    [ST,~]=Data_PadorCut(ST',width_pad,type_pad,flag_pad);
    ST=ST';
    
    %������
    [~,data_mean_ST_ann,data_mean_ST_other,~,~,~]=Data_out(data_p,ST,Cyc,Fre,ind_ann_cyc1);
    
    f_nn=find(FF=='.',1,'last')-1;
    if FlagPanE==1
        %��ͼ������
        data_p(wz_QS)=NaN;
        data_mean_ST_ann(wz_QS)=NaN;
        data_mean_ST_other(wz_QS)=NaN;
        Plot_Ann(time_p,data_p,data_mean_ST_ann,data_mean_ST_other,ann_cyc1,ann_cyc2,PlotLanMode);
        Figname=[[pwd,'\',FigPathU,'\'],FF(1:f_nn),'_�������Ϣ��ȡ'];
        print(gcf,[Figname,'.png'],'-dpng','-r600');
        saveas(gcf,Figname,'pdf');
        close all;
    end
    
    if FlagOanE==1
        %���
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
    
    StaLatLon=[Olon(iiNFZ),Olat(iiNFZ)];%̨վ��γ��
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
        if NaNflag==1%�޵���
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
    YL={'����Ƶ�ι�һ��ƽ������ONA';['(����',num2str(ann_cyc2(2)),'-',num2str(ann_cyc2(1)),'��)']};
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
        if NaNflag==1%�޵���
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
    YL={'��Ƶ�ι�һ��ƽ������ANA';['(����',num2str(ann_cyc1(2)),'-',num2str(ann_cyc1(1)),'��)']};
    OutDZML2=LQPlotSingleCurve1(ZhiBiao2,YL,SelectEC,DisRange,Utime,Udata,outFt);
    LQMolchan1DNew(ZhiBiao2,SelectEC,Utime,Udata,NoThreshold,DisRange,outFt);
    clear SelectEC
    
    outname=strcat(pwd,'\',ResultPathU,'\',FF(1:f_nn),'_BestIndex.mat');
    save(outname,'ZhiBiao1','ZhiBiao2','OutDZML1','OutDZML2');
    clear OutDZML1 OutDZML2
end
%LQMolchanCompute2_ParNew�����˲��������ý���ǰ��ֹMatlab���м��㻷��
delete(gcp('nocreate'));