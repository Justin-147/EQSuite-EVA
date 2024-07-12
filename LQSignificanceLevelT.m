%根据二项分布和正态分布计算显著性水平a
%置信水平=1-a
%作者：刘琦
%单位：中国地震局地震预测研究所
%Compute significance level alpha according to the binomial distribution and normal distribution
%Input: all number of earthquakes,number of correctly predicted earthquakes,space-time occupation
%Output: alpha
%Author: Liu Qi
%Institution: IEF
%Contact: liuqi@ief.ac.cn
%Last Modified: 2022/2/27
function alpha=LQSignificanceLevelT(N,h,tao,mode)
%mode: 1 according to the binomial distribution; 2 according to the binomial distribution and normal distribution
ErrorThresh=1e-5;
if (N==h & tao==1) | (h==0 & tao==0)
    alpha=1;
    return
end
if N>=20 && mode==2%normal distribution
    X=N/(N+1)*(h-N*tao)/sqrt(N*tao*(1-tao));
    alpha1=0;
    alpha2=1;
    alpha3=0.5;
    Xalpha3=icdf('normal',1-alpha3,0,1);
    if X==Inf
        alpha=alpha1;
    elseif X==-Inf
        alpha=alpha2;
    else
        ii=0;
        while abs(X-Xalpha3)>ErrorThresh
            ii=ii+1;
            if Xalpha3>X
                alpha1=alpha3;
            elseif Xalpha3<X
                alpha2=alpha3;
            else
                break;
            end
            alpha3=(alpha1+alpha2)/2;
            Xalpha3=icdf('normal',1-alpha3,0,1);
            if ii>100
                alpha3=NaN;
                break;
            end
        end
        alpha=alpha3;
    end
else
    ii=h:1:N;
    talpha=factorial(N)./factorial(ii)./factorial(N-ii).*tao.^ii.*(1-tao).^(N-ii);
    alpha=sum(talpha);
end
end