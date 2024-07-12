%-------------------------------------------------------------------------
%  dataz=RepInvalidX(dataz,QS,lx)
%  用邻近值简单进行缺数插值
%  作者：刘琦 单位：中国地震局地震预测研究所 最后修改时间2018-10-27 14:06
%-------------------------------------------------------------------------
function [dataz,index2]=RepInvalidX(dataz,QS,lx)
xx=(1:1:length(dataz))';
index2=union(find(dataz==QS),find(isnan(dataz)));
if isempty(index2)%全部不缺数不处理
    return;
end
if length(index2)==length(dataz)%全部为缺数不处理
    return;
end
index=setdiff(xx,index2);
dataz=interp1(index,dataz(index),xx,'nearest');%matlab插值函数要求插值范围不能超过基础数据范围
dataz(index2(index2<index(1)))=dataz(index(1));
dataz(index2(index2>index(end)))=dataz(index(end));
end