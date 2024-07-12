%数据扩边或截边
%2019-02-20，刘琦
function [data_p,oflag]=Data_PadorCut(data_p,width_pad,type_pad,flag_pad)
oflag=1;%顺利完成了操作
switch flag_pad %'pad'扩边,'cut'截边
    case 'pad'
        if ~isempty(data_p)
            switch type_pad
                case 'all'
                    data_p=[repmat(data_p(1,:),width_pad,1);data_p;repmat(data_p(end,:),width_pad,1)];
                case 'left'
                    data_p=[repmat(data_p(1,:),width_pad,1);data_p];
                case 'right'
                    data_p=[data_p;repmat(data_p(end,:),width_pad,1)];
                otherwise
                    oflag=0;
            end
        end
    case 'cut'
        switch type_pad
            case 'all'
                if length(data_p)>2*width_pad
                    data_p=data_p(width_pad+1:end-width_pad,:);
                else
                    oflag=0;
                end
            case 'left'
                if length(data_p)>width_pad
                    data_p=data_p(width_pad+1:end,:);
                else
                    oflag=0;
                end
            case 'right'
                if length(data_p)>width_pad
                    data_p=data_p(1:end-width_pad,:);
                else
                    oflag=0;
                end
            otherwise
                oflag=0;
        end
    otherwise
        oflag=0;
end
end