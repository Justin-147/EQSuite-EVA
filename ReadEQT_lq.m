%�й�����ֵ���Ԥ���о����������ƣ�������ʱ��2022-2-25��liuqi@ief.ac.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ȡ����EQT����Ŀ¼
function [eqDateTime,eqLat,eqLon,eqMag,eqDepth,eqLocation]=ReadEQT_lq(Filename)
eqtStr=textread(Filename, '%s', 'delimiter', '\n', 'whitespace','');
eqtStr=char(eqtStr);
eqDateTime=eqtStr(:,2:15);
eqLat=str2num(eqtStr(:,17:21));
eqLon=str2num(eqtStr(:,23:28));
eqMag=str2num(eqtStr(:,29:31));
eqDepth=str2num(eqtStr(:,34:35));
eqLocation=eqtStr(:,40:end);
end