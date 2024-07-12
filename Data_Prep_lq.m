%数据预处理：补断数、去台阶和突跳、邻近值简单进行缺数插值
%2019-02-20，刘琦
function [time_p,data_p,oflag,wz_QS]=Data_Prep_lq(time_p,data_p,QS)
oflag=1;%输出数据不全为无效数据

[data_p,time_p]=FillGap(data_p,time_p,QS);%填补断数,文件中仍然存在缺数标记
if all(data_p==QS)%全部为缺数
    oflag=0;
    return;
end

lx=3;%1无差别归零,2年归零,3月归零,4日归零
data_p=EraDoubleS(data_p,time_p,QS,lx);%去突跳和台阶

[data_p,wz_QS]=RepInvalidX(data_p,QS,'nearest');%用邻近值简单进行缺数插值,返回了缺数位置
end