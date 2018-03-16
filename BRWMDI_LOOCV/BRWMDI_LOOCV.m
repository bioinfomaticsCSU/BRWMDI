function BRWMDI_LOOCV()
     
%SNF parameters,K:number of neighbors,Nt:Number of Iterations
K = 4;%number of neighbors
Nt = 5; %Number of Iterations
         
%BRWH parameters, alpha:decay facor; 
%Il:maximum iteration numbers of disease network
%Ir:maximum iteration numbers of microbe network
alpha = 0.2;
Il = 1;
Ir = 2;

%predict disease-related microbe based on KATZHMDA in the term of leave-one-out cross validation
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
%implement leave-one-out cross validation
for cv=1:pp 
    %obtain training sample
    load interaction;
    
    interaction(A(cv,1),A(cv,2))=0;
    

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
    
   
    Y = interaction;
    F = BRWMDI(Y, dsymsim, K, Nt, alpha, Il, Ir);

    %obtain the score of tested  disease-microbe interaction
    finalscore=F(A(cv,1),A(cv,2));

    % make the score of seed  disease-microbe interactions as zero
    for i=1:nd
        for j=1:nm
            if interaction(i,j)==1
               F(i,j)=-10000;
            end
        end
    end

    % obtain the position of tested disease-microbe interaction as variable globalposition(1,cv),
    [ll1,mm1]=size(find(F>=finalscore));
    [ll2,mm2]=size(find(F>finalscore));
    globalposition(1,cv)=ll2+1+(ll1-ll2-1)/2;

end
save('globalposition.mat','globalposition');   
end


        
        
        
    
   



