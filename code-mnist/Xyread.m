function [X,y] = Xyread(fea,gnd,type)
[n,~] = size(gnd);
partitional = floor(n/3)/n;
CVO = cvpartition(n,"holdout",partitional);
trvalidx = CVO.training;
testidx = CVO.test;
Xte = fea(testidx,:);
yte = gnd(testidx,:);
Xtrval = fea(trvalidx,:);
ytrval = gnd(trvalidx,:);

[n1,~] = size(ytrval);
partitional2 = 1/2;
CVO2 = cvpartition(n1,"holdout",partitional2);
tridx = CVO2.training;
validx = CVO2.test;
Xtr = Xtrval(tridx,:);
ytr = ytrval(tridx,:);
Xval = Xtrval(validx,:);
yval = ytrval(validx,:);

switch type
    case{'tr'}
        X = Xtr;
        y = ytr;
    case{'val'}
        X = Xval;
        y = yval;
    case{'te'}
        X = Xte;
        y = yte;
end
end






