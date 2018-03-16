function BRWMDI5cv()
    %predict disease-related microbe based on BRWH-MDI in the term of 5-fold cross validation
    %SNF parameters,K:number of neighbors,Nt:Number of Iterations
    K = 4;%number of neighbors
    Nt = 5; %Number of Iterations
         
    %Random walk parameters, alpha:decay facor  0.2
    %Il:maximum iteration numbers of diseases network 
    %Ir:maximum iteration numbers of microbes network 
    alpha = 0.2;
    Il = 1;
    Ir = 2;
    
    %A: Binary relations between disease and microbe, 1st column:disease, 2nd column:microbe
    A=textread('knowndiseasemicrobeinteraction.txt');
    % nd:the number of diseases
    % nm:the number of microbe
    % pp:the number of known diseae-microbe associations
    nd=max(A(:,1)); 
    nm=max(A(:,2));
    [pp,qq]=size(A);
    
    %interaction: adjacency matrix for the disease-microbe association network
    %interaction(i,j)=1 means disease i is related to microbe j
    for i=1:pp
        interaction(A(i,1),A(i,2))=1;
    end
    save interaction interaction;
 
    %obtain the symptom-based similarity of diseases
    load SymtomBased.mat;
    [simnum,colnum]=size(SymtomBased);
    for i=1:simnum
        dsymsim(SymtomBased(i,1),SymtomBased(i,2))=SymtomBased(i,3);
        dsymsim(SymtomBased(i,2),SymtomBased(i,1))=SymtomBased(i,3);
    end
    for i=1:nd
        dsymsim(i,i)=1;
    end    
     
    %implement 5-fold cross validation
    x=randperm(pp)';
    T=1;

    for cv=1:5
        load interaction interaction;
        if cv<5
            B=A(x((cv-1)*floor(pp/5)+1:floor(pp/5)*cv),:);
            %obtain training sample
            for i=1:floor(pp/5)
                interaction(B(i,1),B(i,2))=0;
            end
        else
            B=A(x((cv-1)*floor(pp/5)+1:pp),:);
            %obtain training sample
            for i=1:pp-floor(pp/5)*4
                interaction(B(i,1),B(i,2))=0;
            end
        end
        
        Y = interaction;
        F = BRWMDI(Y, dsymsim, K, Nt, alpha, Il, Ir);
        
        [size1B,size2B]=size(B);
        % obtain the score of tested  disease-microbe interaction
        for i=1:size1B
            finalscore(i,1)=F(B(i,1),B(i,2));
        end
        %make the score of seed  disease-microbe interactions as zero
        for i=1:nd
            for j=1:nm
                if interaction(i,j)==1
                   F(i,j)=-10000;
                end
            end
        end
        for qq=1:size1B
            %obtain the position of tested disease-microbe interaction as variable position(1,cv), 
            [ll1,mm1]=size(find(F>=finalscore(qq)));
            [ll2,mm2]=size(find(F>finalscore(qq)));
            position(1,T)=ll2+1+(ll1-ll2-1)/2;
            T=T+1;
        end

    end
    save('position.mat','position');  

end



        
        
        
    
   



