function z =  accu_20news(y,Xte,yte)
c = 20;
[n,p] = size(Xte);
W = reshape(y,p,c);
Y_pred = Xte*W;
maxY = max((Y_pred'));
y_pred = zeros(n,1);
for i = 1:n
    temp = (Y_pred(i,:)==maxY(i));
    y_pred(i) = find(temp>0);
end
true_pred = (y_pred==yte);
z = sum(true_pred)/n;
end