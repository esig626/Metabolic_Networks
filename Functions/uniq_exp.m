function [expression,expressionid]=uniq_exp(data,ori_id)
%all previous versions are wrong, this is the only corret one. 12-10-2011
expressionid=unique(ori_id);
uniq_id=expressionid;
R=size(uniq_id,1);
expression=zeros(R,size(data,2));
for i=1:R
    ind=ismember(ori_id,uniq_id(i));
    if sum(ind)==1
        expression(i,:)=data(ind,:);
    else
       expression(i,:)=max(data(ind,:));
    end
end

    