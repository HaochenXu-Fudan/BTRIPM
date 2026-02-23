function accuracy = mnist_accuracy(y_bvfim,y_bvfimzj,y_Tregionzj,y_Tregion,input_val,output_val)
[n,m] = size(input_val);
% pressume n is the feature's, m is the the data
[K,~] = size(y_bvfim);
accuracy_bvfim = zeros(K,1);
for k = 1:K
    W = y_bvfim(k,:);
    W_ori = reshape(W(1:end-10),n,10)';
    b = W(end-9:end)';
    y = zeros(m,1);
%     this does not consider the equality exists.
    logit = W_ori*input_val+b;
    mlogit = ones(10,1)*max(logit);
    arg = (logit==mlogit);
    for i = 1:m
        y(i) = (arg(output_val(i)+1,i)==1);
    end
accuracy_bvfim(k) = sum(y)/m;
end

[K,~] = size(y_bvfimzj);
accuracy_bvfimzj = zeros(K,1);
for k = 1:K
    W = y_bvfimzj(k,:);
    W_ori = reshape(W(1:end-10),n,10)';
    b = W(end-9:end)';
    y = zeros(m,1);
%     this does not consider the equality exists.
    logit = W_ori*input_val+b;
    mlogit = ones(10,1)*max(logit);
    arg = (logit==mlogit);
    for i = 1:m
        y(i) = (arg(output_val(i)+1,i)==1);
    end
accuracy_bvfimzj(k) = sum(y)/m;
end

[K,~] = size(y_Tregionzj);
accuracy_Tregionzj = zeros(K,1);
for k = 1:K
    W = y_Tregionzj(k,:);
    W_ori = reshape(W(1:end-10),n,10)';
    b = W(end-9:end)';
    y = zeros(m,1);
%     this does not consider the equality exists.
    logit = W_ori*input_val+b;
    mlogit = ones(10,1)*max(logit);
    arg = (logit==mlogit);
    for i = 1:m
        y(i) = (arg(output_val(i)+1,i)==1);
    end
accuracy_Tregionzj(k) = sum(y)/m;
end

[K,~] = size(y_Tregion);
accuracy_Tregion = zeros(K,1);
for k = 1:K
    W = y_Tregion(k,:);
    W_ori = reshape(W(1:end-10),n,10)';
    b = W(end-9:end)';
    y = zeros(m,1);
%     this does not consider the equality exists.
    logit = W_ori*input_val+b;
    mlogit = ones(10,1)*max(logit);
    arg = (logit==mlogit);
    for i = 1:m
        y(i) = (arg(output_val(i)+1,i)==1);
    end
accuracy_Tregion(k) = sum(y)/m;
end


accuracy = struct();
accuracy.bvfim = accuracy_bvfim;
accuracy.bvfimzj = accuracy_bvfimzj;
accuracy.Tregionzj = accuracy_Tregionzj;
accuracy.Tregion = accuracy_Tregion;
end