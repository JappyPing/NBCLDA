function NBCLDA_LOOCV(A,B,D,E,F,G,C,H)
%%
%----------------------------------------------------------------------------
% Note: The implementation of this function refers to the code in Chen et al.'s 
% work(DOI: https://academic.oup.com/bioinformatics/article/33/5/733/2736222)  
%----------------------------------------------------------------------------
%%
% A: adjacency matrix for the miRNA-lncRNA associations
% B: adjacency matrix for the miRNA-disease associations
% D: adjacency matrix for the gene-lncRNA associations
% E: adjacency matrix for the gene-disease associations
% F: adjacency matrix for the miRNA-gene associations
% C: adjacency matrix for the lncRNA-disease associations
% ld records the location of non-zero elements in LD
% DS represents the disease semantic similarity

interaction=C;
% pp:the number of known lncRNA-diseae associations
[pp,qq]=size(G);
save interaction interaction;
% implement leave-one-out cross validation
for cv=1:pp 
    cv
% obtain training sample
    load interaction;
    M=interaction;
    M(G(cv,1),G(cv,2))=0;
    [score] = NBCLDA_GN2_SD(A,B,M,D,E,F,H);
% obtain the score of tested LncRNA-Disease interaction
    finalscore=score(G(cv,2),G(cv,1));
% make the score of seed  LncRNA-Disease interactions as zero
    [nd,nm]=size(M);
for i=1:nd
    for j=1:nm
        if interaction(i,j)==1
           score(j,i)=-10000;
        end
    end
end
% obtain the position of tested LncRNA-Disease interaction as variable globalposition(1,cv),
[ll1,mm1]=size(find(score>=finalscore));
[ll2,mm2]=size(find(score>finalscore));
globalposition(1,cv)=ll2+1+(ll1-ll2-1)/2;
end
save('globalposition.mat','globalposition');   
end


        
        
        
    
   



