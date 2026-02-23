function [input_val, output_val] = dataforval(funname,m)
switch funname
    case{'mnist5000'}
        mnist_image = loadMNISTImages('data/mnist_train_image.idx3-ubyte');
        mnist_label = loadMNISTLabels('data/mnist_train_label.idx1-ubyte');
        R = randperm(60000);
        input_val = mnist_image(:,R(1:m));
        output_val = mnist_label(R(1:m),1);
    case{'fashion_mnist'}
        fmnist_image = loadMNISTImages('data/fashion_train-images-idx3-ubyte');
        fmnist_label = loadMNISTLabels('data/fashion_train-labels-idx1-ubyte');
        R = randperm(60000);
        input_val = fmnist_image(:,R(1:m));
        output_val = fmnist_label(R(1:m),1);
end
end