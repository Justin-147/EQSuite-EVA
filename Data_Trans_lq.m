%�������ݲ�ת���ɱ�׼2�и�ʽ
%2019-02-19������
%��Ҫ�㷨��������д��PAS_FD_Prep_Trans_Callback����
function [timet,dataz]=Data_Trans_lq(dbfile)
tmp=load(dbfile); [~,N]=size(tmp);
if N==2%��׼��ʽ
    timet=tmp(:,1);
    dataz=tmp(:,2);
elseif N==25%����ֵ
    dataz=tmp(:,2:25);    timet=tmp(:,1);
    dataz=dataz';
    dataz=dataz(:);
    aa=timet*100;
    aa=repmat(aa,1,24);
    bb=0:1:23;
    bb=repmat(bb,size(aa,1),1);
    timet=aa+bb;
    timet=timet';
    timet=timet(:);
elseif N==1441%����ֵ
    dataz=tmp(:,2:1441);    timet=tmp(:,1);
    dataz=dataz';
    dataz=dataz(:);
    aa=timet*10000;
    aa=repmat(aa,1,1440);
    bb=0:100:2300;
    bb=repmat(bb,60,1);
    bb=repmat(bb(:)',size(aa,1),1);
    cc=0:1:59;
    cc=repmat(cc,size(aa,1),24);
    timet=aa+bb+cc;
    timet=timet';
    timet=timet(:);
else
end
end