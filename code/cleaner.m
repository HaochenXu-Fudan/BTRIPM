function FI = cleaner(c,x_bvfim,x_bvfimzj,x_Tregionzj,x_Tregion,thr)
[K,~] = size(x_bvfim);
FI_bvfim = zeros(K,1);
boundline = thr;
for i = 1:K
    x = x_bvfim(i,:);
    tp = sum(x(c)<boundline);
    fp = sum(x<boundline)-tp;
    fn = sum(x(c)>=boundline);
    tn = sum(x>=boundline)-fn;
    prec = tp/(tp+fp);
    rec = tp/(tp+fn);
    FI_bvfim(i) = 2*prec*rec/(prec+rec);
end
[K,~] = size(x_bvfimzj);
FI_bvfimzj = zeros(K,1);
% boundline = thr.pro;
for i = 1:K
    x = x_bvfimzj(i,:);
    tp = sum(x(c)<boundline);
    fp = sum(x<boundline)-tp;
    fn = sum(x(c)>=boundline);
    tn = sum(x>=boundline)-fn;
    prec = tp/(tp+fp);
    rec = tp/(tp+fn);
    FI_bvfimzj(i) = 2*prec*rec/(prec+rec);
end
[K,~] = size(x_Tregionzj);
FI_Tregionzj = zeros(K,1);
boundline = thr;
for i = 1:K
    x = x_Tregionzj(i,:);
    tp = sum(x(c)<boundline);
    fp = sum(x<boundline)-tp;
    fn = sum(x(c)>=boundline);
    tn = sum(x>=boundline)-fn;
    prec = tp/(tp+fp);
    rec = tp/(tp+fn);
    FI_Tregionzj(i) = 2*prec*rec/(prec+rec);
end

[K,~] = size(x_Tregion);
FI_Tregion = zeros(K,1);
% boundline = thr.Tregion;
for i = 1:K
    x = x_Tregion(i,:);
    tp = sum(x(c)<boundline);
    fp = sum(x<boundline)-tp;
    fn = sum(x(c)>=boundline);
    tn = sum(x>=boundline)-fn;
    prec = tp/(tp+fp);
    rec = tp/(tp+fn);
    FI_Tregion(i) = 2*prec*rec/(prec+rec);
end
FI = struct();
FI.bvfim = FI_bvfim;
FI.bvfimzj = FI_bvfimzj;
FI.Tregionzj = FI_Tregionzj;
FI.Tregion = FI_Tregion;
end