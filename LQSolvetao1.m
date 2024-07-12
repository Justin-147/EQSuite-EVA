%中国地震局地震预测研究所刘琦编制，最后调试时间2022-3-3，liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%给定显著性水平求解对应时空占有率tao
%tao异常的时空占有率，发出预测“警报”的时空范围与总的时空范围之比
%N总地震数
%h预测“有震”而实际发震的地震数，击中数
%a显著性水平
function tao=LQSolvetao1(N,h,a)
%nchoosek可以直接计算组合数，但考虑不用循环，放弃使用
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
    if a3*a1>0&&a3*a2<0%与a1同号，与a2异号
        tao1=tao3;
        tao3=(tao2+tao1)/2;
        a1=LQSignificanceLevelT(N,h,tao1,2)-a;
        a3=LQSignificanceLevelT(N,h,tao3,2)-a;
    elseif a3*a1<0&&a3*a2>0%与a2同号，与a1异号
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