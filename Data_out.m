%�����һ������
%2022-02-25������
function [data_mean_ST,data_mean_ST_ann,data_mean_ST_other,Fre_max_ST,Fre_max_ST_ann,Fre_max_ST_other]=Data_out(data_p,ST,Cyc,Fre,ind_ann_cyc1)
del_factor=0.05;
del_len=round(length(data_p)*del_factor);%ȥ��һ�������������Сֵ��ʹ��ƽ��ֵ������Ը�����
%%%%%%%%%
data_mean_ST=mean(abs(ST));%����Ƶ��
data_mean_ST=Data_Norm(data_mean_ST,del_len);%��ֵ��һ��
%
data_mean_ST_ann=mean(abs(ST(ind_ann_cyc1(1):ind_ann_cyc1(2),:)));%������Ƶ��
data_mean_ST_ann=Data_Norm(data_mean_ST_ann,del_len);%��ֵ��һ��
%
indt=1:length(Cyc);
indt=setdiff(indt,ind_ann_cyc1(1):ind_ann_cyc1(2));
data_mean_ST_other=mean(abs(ST(indt,:)));%��������Ƶ��
data_mean_ST_other=Data_Norm(data_mean_ST_other,del_len);%��ֵ��һ��
%%%%%%%%%
[~,ind]=max(abs(ST));%����Ƶ��
Fre_max_ST=Fre(ind);
Fre_max_ST=Data_Norm(Fre_max_ST,del_len);%�����ȶ�Ӧ��Ƶ�����й�һ��
%
[~,ind]=max(abs(ST(ind_ann_cyc1(1):ind_ann_cyc1(2),:)));%������Ƶ��
Fre_max_ST_ann=Fre(ind);
Fre_max_ST_ann=Data_Norm(Fre_max_ST_ann,del_len);%�����ȶ�Ӧ��Ƶ�����й�һ��
%
[~,ind]=max(abs(ST(indt,:)));%��������Ƶ��
Fre_max_ST_other=Fre(ind);
Fre_max_ST_other=Data_Norm(Fre_max_ST_other,del_len);%�����ȶ�Ӧ��Ƶ�����й�һ��
end

%ȥ��ֵ����һ��
function data_out=Data_Norm(data_input,del_len)
tmp=sort(data_input);
tmp=mean(tmp(del_len+1:end-del_len));
data_out=(data_input-tmp)/tmp;
end