%����Ԥ������������ȥ̨�׺�ͻ�����ڽ�ֵ�򵥽���ȱ����ֵ
%2019-02-20������
function [time_p,data_p,oflag,wz_QS]=Data_Prep_lq(time_p,data_p,QS)
oflag=1;%������ݲ�ȫΪ��Ч����

[data_p,time_p]=FillGap(data_p,time_p,QS);%�����,�ļ�����Ȼ����ȱ�����
if all(data_p==QS)%ȫ��Ϊȱ��
    oflag=0;
    return;
end

lx=3;%1�޲�����,2�����,3�¹���,4�չ���
data_p=EraDoubleS(data_p,time_p,QS,lx);%ȥͻ����̨��

[data_p,wz_QS]=RepInvalidX(data_p,QS,'nearest');%���ڽ�ֵ�򵥽���ȱ����ֵ,������ȱ��λ��
end