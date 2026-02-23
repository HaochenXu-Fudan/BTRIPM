%% IMBLO
% This is a file to demo various interior methods on bilevel optimization
% problems.
%% Function Initialization
clear;
experiment = 5; 
switch experiment
    case 1
        funname = 'toy0';x0 = 3;y0 = 3;x0str = 'x_3';y0str = 'y_3';
    case 2
        funname = 'wifi';x0 = 1.5;y0 = -.1*ones(7,1);
    case 3
        funname = 'crowdsourced';x0 = 2;y0 = ones(28,1);
    case 4
        funname = 'mnist5000';x0 = zeros(5000,1);y0 = .1*ones(7850,1);
    case 5
        funname = 'fashion_mnist';x0 = zeros(5000,1);y0 = .1*ones(7850,1);
end
[F,f,diff,others] = fun_read(funname);
%% bvfim
opts = Opts(funname,'bvfim');
[x_bvfim,y_bvfim,hist_bvfim,time_bvfim] = bvfim(x0,y0,F,f,diff,opts);
%% pro (ABANDONED)
% opts = Opts(funname,'pro');
% [x_pro,y_pro,hist_pro,time_pro] = pro(x0,y0,F,f,diff,opts);
%% rhg & trhg
if experiment<4
    opts = Opts(funname,'rhg');
    [x_rhg,y_rhg,hist_rhg,time_rhg] = rhgtest(x0,y0,F,f,diff,opts);
    if experiment
        opts = Opts(funname,'trhg');
        [x_trhg,y_trhg,hist_trhg,time_trhg] = rhgtest(x0,y0,F,f,diff,opts);
    end
else
    opts = Opts(funname,'bvfimzj');
    [x_bvfimzj,y_bvfimzj,hist_bvfimzj,time_bvfimzj] = bvfim(x0,y0,F,f,diff,opts);
end
%% cg
if experiment<4
    opts = Opts(funname,'cg');
    [x_cg,y_cg,hist_cg,time_cg] = cgtest(x0,y0,F,f,diff,opts);
else
    opts = Opts(funname,'Tregionzj');
    [x_Tregionzj,y_Tregionzj,hist_Tregionzj,time_Tregionzj] = Tregion(x0,y0,F,f,diff,opts);
end
%% Trust region
opts = Opts(funname,'Tregion');
[x_Tregion,y_Tregion,hist_Tregion,time_Tregion] = Tregion(x0,y0,F,f,diff,opts);
%% error compute
% Only in toy example, error = |value-real_val|, in real world datasets, we
% directly use value.
plot_file = 1;
if strcmp(funname,'toy0')||strcmp(funname,'toy2')
    error_bvfim = errorcompute(funname,hist_bvfim);
    error_rhg = errorcompute(funname,hist_rhg);
%     error_trhg = errorcompute(funname,hist_trhg);
    error_cg = errorcompute(funname,hist_cg);
    error_Tregion = errorcompute(funname,hist_Tregion);
    plot_file = 0;
end
%% save result
saldir = './result/';
savePath = [saldir funname '.mat'];
if experiment<4
    if experiment==1
        save(savePath,'hist_Tregion','hist_bvfim','hist_rhg','hist_cg','time_Tregion','time_bvfim','time_rhg','time_cg');        
    else
        save(savePath,'hist_Tregion','hist_bvfim','hist_rhg','hist_trhg','hist_cg','time_Tregion','time_bvfim','time_rhg','time_trhg','time_cg');       
    end
else
    save(savePath,'hist_Tregion','hist_bvfim','hist_Tregionzj','hist_bvfimzj','time_Tregion','time_bvfim','time_Tregionzj','time_bvfimzj');
end
%% plot
plotyes = 1;
if plotyes==1
    if plot_file
        F0 = F(x0,y0);
        f0 = f(x0,y0);
        if experiment<4
            plot_val(funname,hist_bvfim,hist_Tregion,hist_rhg,hist_trhg,hist_cg,F0,f0)
        else
            plot_valmnist(funname,hist_bvfim,hist_Tregion,hist_bvfimzj,hist_Tregionzj,F0,f0)
        end
    else
        plot_error(error_bvfim,error_rhg,error_cg,error_Tregion,funname,x0str,y0str)
        savePath = [saldir funname '_error' '.mat'];
        save(savePath,'error_Tregion','error_bvfim','error_rhg','error_cg');
    end
end
%% Accuracy&FI
plotyes = 1;
if  strcmp(funname,'mnist5000')||strcmp(funname,'fashion_mnist')
    [input_test,output_test] = dataforval(funname,10000);
    accuracy = mnist_accuracy(hist_bvfim.y,hist_bvfimzj.y,hist_Tregionzj.y,hist_Tregion.y,input_test,output_test);
    thr = 0;
    FI = cleaner(others.c,hist_bvfim.x,hist_bvfimzj.x,hist_Tregionzj.x,hist_Tregion.x,thr);
    savePath = [saldir funname '_accu' '.mat'];
    save(savePath,'accuracy','FI');
    if plotyes
        FI.plot_yes = 1;
        plot_accuracy(funname,accuracy,FI);
    end
end