clear;
% filename = '/Users/severe/Downloads/high/20Newsgroups.mat';
load('20Newsgroups.mat');
%% generate [X,y]
[Xtr,ytr] = Xyread(fea,gnd,'tr');
[Xval,yval] = Xyread(fea,gnd,'val');
[Xte,yte] = Xyread(fea,gnd,'te');
save('Xtr','Xtr');
save('Xval','Xval');
save('Xte','Xte');
save('yte','yte');
n = length(ytr); c = 20;
Ytr = zeros(n,c);
for i = 1: n
    Ytr(i,ytr(i)) = 1;
end
n = length(yval); c = 20;
Ytr = sparse(Ytr);
Yval = zeros(n,c);
for i = 1: n
    Yval(i,yval(i)) = 1;
end
Yval = sparse(Yval);
save('Ytr','Ytr');
save('Yval','Yval');