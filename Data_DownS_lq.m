%数据降采样
%2019-02-20，刘琦
%主要算法据刘琦编写的PAS_FD_Prep_DownS_Callback函数
function [time_p,data_p,oflag]=Data_DownS_lq(time_p,data_p,lx)
oflag=1;%完成了数据降采样
ODL=length(num2str(time_p(1)));
switch lx
    case '整日值'
        if ODL==8%整日值
            return;
        elseif ODL==10%整时值
            CXS=100;
        elseif ODL==12%分钟值
            CXS=10000;
        elseif ODL==14%秒值
            CXS=1000000;
        else
            oflag=0;
            return;
        end
        ind=mod(time_p,CXS)==0;
        data_p=data_p(ind);
        time_p=time_p(ind)/CXS;        
    case '整时值'
        if ODL==10%整时值
            return;
        elseif ODL==12%分钟值
            CXS=100;
        elseif ODL==14%秒值
            CXS=10000;
        else
            oflag=0;
            return;
        end
        ind=mod(time_p,CXS)==0;
        data_p=data_p(ind);
        time_p=time_p(ind)/CXS;
    case '分钟值'
        if ODL==12%分钟值
            return;
        elseif ODL==14%秒值
            CXS=100;
        else
            oflag=0;
            return;
        end
        ind=mod(time_p,CXS)==0;
        data_p=data_p(ind);
        time_p=time_p(ind)/CXS;
    otherwise
        oflag=0;
end
end