function result = Measures(Y, predY)
% result = [ACC, NMI, ARI, F-score];
if size(Y,2) ~= 1
    Y = Y';
end;
if size(predY,2) ~= 1
    predY = predY';
end;

% bestMap
predY = bestMap(Y, predY);
if size(Y)~=size(predY)
    predY=predY';
end
error_cnt = sum(Y ~= predY);
AC = length(find(Y == predY))/length(Y);
[~,nmi_value,~] = compute_nmi(Y', predY');
[ARI,~,~,~] = valid_RandIndex(Y', predY');
f = compute_f(Y, predY);
result = [AC, nmi_value, ARI, f,error_cnt];