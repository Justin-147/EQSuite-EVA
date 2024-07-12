%根据二项分布和正态分布计算R0
%矩阵思想
%作者：刘琦
%单位：中国地震局地震预测研究所
%Compute R0 according to the binomial distribution and normal distribution
%Input: successAlamCount, missAlamCount, [alpha], [mode]
%Output: R0
%Author: Liu Qi
%Institution: IEF
%Contact: liuqi@ief.ac.cn
%Last Modified: 2022/2/22
function [R0]=LQGetR0(varargin)
R0=NaN;
P=0:0.00001:1;
alpha=0.025;%97.25% confidence interval
mode=2;%1 according to the binomial distribution; 2 according to the binomial distribution and normal distribution
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
        X=n/(n+1)*(k-n*P)./sqrt(n*P.*(1-P));
        lessAlpha=abs(X-Xalpha);%X≥Xalpha
        [~,II]=min(lessAlpha);
        R0=k/n-P(II);
    else
        P=repmat(P,n-k+1,1);
        ii=k:1:n;
        ii=repmat(ii',1,size(P,2));
        talpha=factorial(n)./factorial(ii)./factorial(n-ii).*P.^ii.*(1-P).^(n-ii);
        Talpha=sum(talpha,1);
        lessAlpha=abs(Talpha-alpha);%≤alpha
        [~,II]=min(lessAlpha);
        R0=k/n-P(1,II);
    end
end
end