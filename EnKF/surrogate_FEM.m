function [Imax,Imin,mse,b]=surrogate_FEM(Nc0,EIs)
%% SURROGATE TO ESTIMATE THE ERROR INDICATOR
    
    x1 = Nc0(:,1);
    x2 = Nc0(:,2);
    x3 = Nc0(:,3);

    y = EIs';
    
    X = [ones(size(x1)), x1, x2, x3, x1.*x2, x1.*x3, x2.*x3, (x1.*x2).*x3, ...
        x1.^2, x2.^2, x3.^2];
%     X = [ones(size(x1)), x1, x2, x3, x1.*x2, x1.*x3, x2.*x3, (x1.*x2).*x3];
    [b,~,r] = regress(y,X);
    mse = (sum(r.^2))/size(r,1);
    
    x1fit = linspace(50e6,50e8,10);
    x2fit = linspace(0.02,0.11,10); 
    x3fit = linspace(0.002,0.011,10);
    [X1FIT,X2FIT,X3FIT] = meshgrid(x1fit,x2fit,x3fit);
    
    YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X3FIT + b(5)*X1FIT.*X2FIT + ...
        b(6)*X1FIT.*X3FIT + b(7)*X2FIT.*X3FIT + b(8)*(X1FIT.*X2FIT).*X3FIT + ...
        (b(9)*(X1FIT.^2)) + (b(10)*(X2FIT.^2)) + (b(11)*(X3FIT.^2));
    
    
%     YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X3FIT + b(5)*X1FIT.*X2FIT + ...
%         b(6)*X1FIT.*X3FIT + b(7)*X2FIT.*X3FIT + b(8)*(X1FIT.*X2FIT).*X3FIT;
    Imax = max(YFIT,[],'all','linear');
    Imin = min(YFIT,[],'all','linear');    

end

