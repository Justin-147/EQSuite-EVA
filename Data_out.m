%结果归一化整理
%2022-02-25，刘琦
function [data_mean_ST,data_mean_ST_ann,data_mean_ST_other,Fre_max_ST,Fre_max_ST_ann,Fre_max_ST_other]=Data_out(data_p,ST,Cyc,Fre,ind_ann_cyc1)
del_factor=0.05;
del_len=round(length(data_p)*del_factor);%去掉一定比例的最大、最小值后，使得平均值计算相对更合理
%%%%%%%%%
data_mean_ST=mean(abs(ST));%整个频带
data_mean_ST=Data_Norm(data_mean_ST,del_len);%幅值归一化
%
data_mean_ST_ann=mean(abs(ST(ind_ann_cyc1(1):ind_ann_cyc1(2),:)));%年周期频带
data_mean_ST_ann=Data_Norm(data_mean_ST_ann,del_len);%幅值归一化
%
indt=1:length(Cyc);
indt=setdiff(indt,ind_ann_cyc1(1):ind_ann_cyc1(2));
data_mean_ST_other=mean(abs(ST(indt,:)));%其它周期频带
data_mean_ST_other=Data_Norm(data_mean_ST_other,del_len);%幅值归一化
%%%%%%%%%
[~,ind]=max(abs(ST));%整个频带
Fre_max_ST=Fre(ind);
Fre_max_ST=Data_Norm(Fre_max_ST,del_len);%最大幅度对应的频率序列归一化
%
[~,ind]=max(abs(ST(ind_ann_cyc1(1):ind_ann_cyc1(2),:)));%年周期频带
Fre_max_ST_ann=Fre(ind);
Fre_max_ST_ann=Data_Norm(Fre_max_ST_ann,del_len);%最大幅度对应的频率序列归一化
%
[~,ind]=max(abs(ST(indt,:)));%其它周期频带
Fre_max_ST_other=Fre(ind);
Fre_max_ST_other=Data_Norm(Fre_max_ST_other,del_len);%最大幅度对应的频率序列归一化
end

%去均值并归一化
function data_out=Data_Norm(data_input,del_len)
tmp=sort(data_input);
tmp=mean(tmp(del_len+1:end-del_len));
data_out=(data_input-tmp)/tmp;
end