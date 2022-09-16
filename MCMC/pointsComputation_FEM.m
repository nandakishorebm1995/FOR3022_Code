function Ncn = pointsComputation_FEM(b,Imax,Imin,mse,params)
%% Identifies parameter sample points with maximum error indicator values 
%%
    nT = 500;                                             
    trialPoints = params(randperm(size(params,1),nT),:);
    nTau = linspace(0.1,3,10);                            
    mu_star = zeros(length(nTau),3);
    
    for j = 1:length(nTau)
        
        tau = Imax+(nTau(j)*(Imax-Imin));        
        yfit = b(1) + b(2).*trialPoints(:,1) + b(3).*trialPoints(:,2) +...
               b(4).*trialPoints(:,3) + (b(5).*trialPoints(:,1)).*trialPoints(:,2) + ...
               (b(6).*trialPoints(:,1)).*trialPoints(:,3) + (b(7).*trialPoints(:,2)).*trialPoints(:,3) +...
               ((b(8).*trialPoints(:,1)).*trialPoints(:,2)).*trialPoints(:,3) + ...
               (b(9).*(trialPoints(:,1).^2)) + (b(10).*(trialPoints(:,2).^2)) + ...
               (b(11).*(trialPoints(:,1).^2));   
           
%         yfit = b(1) + b(2).*trialPoints(:,1) + b(3).*trialPoints(:,2) + b(4).*trialPoints(:,3)+ (b(5).*trialPoints(:,1)).*trialPoints(:,2) + ...
%                (b(6).*trialPoints(:,1)).*trialPoints(:,3) + (b(7).*trialPoints(:,2)).*trialPoints(:,3) +...
%                ((b(8).*trialPoints(:,1)).*trialPoints(:,2)).*trialPoints(:,3) ;
        ctaus=zeros(1,nT);
       
       for k = 1:nT             
           ptau=(yfit(k)-tau)/sqrt(mse);
           ctau=normpdf(ptau,yfit(k),sqrt(mse)); 
           ctaus(k)=ctau;          
       end
       
       [~,idxctaus] = max(ctaus);
       mu_star(j,:) = trialPoints(idxctaus,:); 
       
    end
        
    myfunc = @(XX,KK)(kmeans(XX, KK, 'emptyaction','singleton','replicate',5));   
    eva = evalclusters(mu_star,myfunc,'CalinskiHarabasz','klist',[1:20]);
    OptClust = eva.OptimalK;

    if isnan(OptClust)
      OptClust = 7;
    end
    
    fileID = fopen('Adaptive_Log_File.txt','a+');
    fprintf(fileID,'Cs: %s \n',string(OptClust));    
    fclose(fileID);
    
    [~,Ncn,~] = kmeans(mu_star,OptClust);
    Ncn(:,2:3) = Ncn(:,2:3);
    
end

