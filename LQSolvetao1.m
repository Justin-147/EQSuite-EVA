%�й�����ֵ���Ԥ���о����������ƣ�������ʱ��2022-3-3��liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����������ˮƽ����Ӧʱ��ռ����tao
%tao�쳣��ʱ��ռ���ʣ�����Ԥ�⡰��������ʱ�շ�Χ���ܵ�ʱ�շ�Χ֮��
%N�ܵ�����
%hԤ�⡰���𡱶�ʵ�ʷ���ĵ�������������
%a������ˮƽ
function tao=LQSolvetao1(N,h,a)
%nchoosek����ֱ�Ӽ���������������ǲ���ѭ��������ʹ��
% if h==0
%     tao=0;
%     return;
% elseif h==N
%     tao=1;
%     return;
% else
%     %errordlg('Can not solve!!');
% end

tao1=0;tao2=1;tao3=(tao2+tao1)/2;
delta=1E-5;

a1=LQSignificanceLevelT(N,h,tao1,2)-a;
a2=LQSignificanceLevelT(N,h,tao2,2)-a;
a3=LQSignificanceLevelT(N,h,tao3,2)-a;

if abs(a1)<=delta
    tao=tao1;
    return;
elseif abs(a2)<=delta
    tao=tao2;
    return;
elseif abs(a3)<=delta
    tao=tao3;
    return;
else
end

while(abs(a1-a2)>delta)
    if a3*a1>0&&a3*a2<0%��a1ͬ�ţ���a2���
        tao1=tao3;
        tao3=(tao2+tao1)/2;
        a1=LQSignificanceLevelT(N,h,tao1,2)-a;
        a3=LQSignificanceLevelT(N,h,tao3,2)-a;
    elseif a3*a1<0&&a3*a2>0%��a2ͬ�ţ���a1���
        tao2=tao3;
        tao3=(tao2+tao1)/2;
        a2=LQSignificanceLevelT(N,h,tao2,2)-a;
        a3=LQSignificanceLevelT(N,h,tao3,2)-a;
    else
        %errordlg('Can not solve!');
        tao=NaN;
        return;
    end
    if abs(a1)<=delta
        tao=tao1;
        return;
    elseif abs(a2)<=delta
        tao=tao2;
        return;
    elseif abs(a3)<=delta
        tao=tao3;
        return;
    else
        continue;
    end
end
end