%�й�����ֵ���Ԥ���о����������ƣ�������ʱ��2022-2-25��liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������оࡢ�𼶡�ʱ��ɸѡ����Ŀ¼
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ind,UseddzDate,UseddzDateNum,UseddzLatitude,UseddzLongitude,UseddzMagnitude,UseddzDepth,UseddzLocation]=LQdzSelect1(LatStation,LonStation,dzDate,dzLatitude,dzLongitude,dzMagnitude,dzDepth,dzLocation,DisRange,MM,DateRangeNum)
ind=true(length(dzLatitude),1);
lendzDate=length(num2str(dzDate(1)));
dzDateNum=datenum(num2str(floor(dzDate/(10^(lendzDate-8)))),'yyyymmdd');
if ~isempty(DisRange)
    dist=vdist(LatStation*ones(length(dzLatitude),1),LonStation*ones(length(dzLatitude),1),dzLatitude,dzLongitude)/1000;%���о࣬��λkm
    ind=ind&dist>=DisRange(1)&dist<DisRange(2);
end
if ~isempty(MM)
    ind=ind&dzMagnitude>=MM(1)&dzMagnitude<MM(2);
end
if ~isempty(DateRangeNum)
    ind=ind&dzDateNum>DateRangeNum(1)&dzDateNum<=DateRangeNum(2);
end
UseddzDate=dzDate(ind);
UseddzDateNum=dzDateNum(ind);
UseddzLatitude=dzLatitude(ind);
UseddzLongitude=dzLongitude(ind);
UseddzMagnitude=dzMagnitude(ind);
UseddzDepth=dzDepth(ind);
UseddzLocation=dzLocation(ind,:);
end