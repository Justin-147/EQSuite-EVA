%根据二项分布和正态分布计算R0
%二分法思想
%作者：刘琦
%单位：中国地震局地震预测研究所
%Compute R0 according to the binomial distribution and normal distribution
%Input: successAlamCount, missAlamCount, [alpha], [mode]
%Output: R0
%Author: Liu Qi
%Institution: IEF
%Contact: liuqi@ief.ac.cn
%Last Modified: 2022/2/27
function [R0]=LQGetR00(varargin)
R0=NaN;
alpha=0.025;%97.25% confidence interval
mode=2;%1 according to the binomial distribution; 2 according to the binomial distribution and normal distribution
ErrorThresh=1e-5;
if nargin<2
    return;%the number of input is insufficient
else
    if nargin>=3
        alpha=varargin{3};%confidence interval
    end
    if nargin>=4
        mode=varargin{4};%computation mode
    end
    k=varargin{1};%successAlamCount
    n=k+varargin{2};%successAlamCount+missAlamCount
    if k==0%successAlamCount==0
        R0=0;
        return;
    end
    if n>=20 & mode==2%normal distribution
        Xalpha=icdf('normal',1-alpha,0,1);
        P1=0;
        P2=1;
        P3=0.5;
        Xalpha3=n/(n+1)*(k-n*P3)/sqrt(n*P3*(1-P3));
        if Xalpha==Inf
            R0=k/n-P1;
        elseif Xalpha==-Inf
            R0=k/n-P2;
        else
            mk=0;
            while abs(Xalpha-Xalpha3)>ErrorThresh
                mk=mk+1;
                if Xalpha3>Xalpha
                    P1=P3;
                elseif Xalpha3<Xalpha
                    P2=P3;
                else
                    break;
                end
                P3=(P1+P2)/2;
                Xalpha3=n/(n+1)*(k-n*P3)/sqrt(n*P3*(1-P3));
                if mk>100
                    P3=NaN;
                    break;
                end
            end
            R0=k/n-P3;
        end
    else
        ii=k:1:n;
        P1=0;
        P2=1;
        P3=0.5;
        if alpha==0
            R0=k/n-P1;
        elseif alpha==1
            R0=k/n-P2;
        else
            alpha3=sum(factorial(n)./factorial(ii)./factorial(n-ii).*P3.^ii.*(1-P3).^(n-ii));
            mk=1;
            while abs(alpha-alpha3)>ErrorThresh
                mk=mk+1;
                if alpha3>alpha
                    P2=P3;
                elseif alpha3<alpha
                    P1=P3;
                else
                    break;
                end
                P3=(P1+P2)/2;
                alpha3=sum(factorial(n)./factorial(ii)./factorial(n-ii).*P3.^ii.*(1-P3).^(n-ii));
                if mk>100
                    P3=NaN;
                    break;
                end
            end
            R0=k/n-P3;
        end
    end
end
end