function [Y_pred,X_pred_idx]=interpLinear(X,Y,X_pred)
% a simple 1D linear interpolation
% Y can be m x n matrix, n is variable
%
[X,idx]=sort(X);
Y=Y(idx,:);
Y_pred=zeros(length(X_pred),size(Y,2));
X_pred_idx=zeros(length(X_pred),1);
for x_idx=1:length(X_pred)
    x_pred=X_pred(x_idx);
    num=length(X);
    idx=num; % search start from last one, find out X samll than x
    while ((idx > 1) && (X(idx) > x_pred))
        idx=idx-1;
    end

    if (idx == num)
        % out interp
        Y_pred(x_idx,:)=(Y(num,:)-Y(num-1,:))/(X(num)-X(num-1))*(x_pred-X(num))+Y(num,:);
        X_pred_idx(x_idx)=idx;
    elseif (idx == 0)
        Y_pred(x_idx,:)=(Y(2,:)-Y(1,:))/(X(2)-X(1))*(x_pred-X(1))+Y(1,:);
        X_pred_idx(x_idx)=idx;
    else
        % linear interpolation
        Y_pred(x_idx,:)=Y(idx,:)+...
            (Y(idx+1,:)-Y(idx,:))*...
            (x_pred-X(idx,:))/...
            (X(idx+1)-X(idx));
        X_pred_idx(x_idx)=idx;
    end
end

end
