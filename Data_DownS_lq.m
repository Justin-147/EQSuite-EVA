%���ݽ�����
%2019-02-20������
%��Ҫ�㷨��������д��PAS_FD_Prep_DownS_Callback����
function [time_p,data_p,oflag]=Data_DownS_lq(time_p,data_p,lx)
oflag=1;%��������ݽ�����
ODL=length(num2str(time_p(1)));
switch lx
    case '����ֵ'
        if ODL==8%����ֵ
            return;
        elseif ODL==10%��ʱֵ
            CXS=100;
        elseif ODL==12%����ֵ
            CXS=10000;
        elseif ODL==14%��ֵ
            CXS=1000000;
        else
            oflag=0;
            return;
        end
        ind=mod(time_p,CXS)==0;
        data_p=data_p(ind);
        time_p=time_p(ind)/CXS;        
    case '��ʱֵ'
        if ODL==10%��ʱֵ
            return;
        elseif ODL==12%����ֵ
            CXS=100;
        elseif ODL==14%��ֵ
            CXS=10000;
        else
            oflag=0;
            return;
        end
        ind=mod(time_p,CXS)==0;
        data_p=data_p(ind);
        time_p=time_p(ind)/CXS;
    case '����ֵ'
        if ODL==12%����ֵ
            return;
        elseif ODL==14%��ֵ
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