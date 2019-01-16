function overallauc = positiontooverallauc(LD,ld)
%%
%----------------------------------------------------------------------------
% Note: The implementation of this function refers to the code in Chen et al.'s 
% work (DOI: https://academic.oup.com/bioinformatics/article/33/5/733/2736222).
%----------------------------------------------------------------------------
%%
% A: adjacency matrix for the lncRNA-disease associations
% B records the location of non-zero elements in LD

load globalposition.mat;
[n,m]=size(A);
[pp,qq]=size(B);
for i=1:pp
    if globalposition(i)>m*n-pp+1
       globalposition(i)=m*n-pp+1;
    end
end
for k=1:m*n-pp+1 
    tp=0;
    for t=1:pp
        if globalposition(1,t)<=k
           tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp;
    fp=k*pp-tp;
    fpr(1,k)=fp/(pp*(m*n-pp));   
end           
plot(fpr,tpr)
hold on
clear area;
area(1,1)=tpr(1,1)*fpr(1,1)/2;
for k=2:m*n-pp+1
    area(1,k)=[tpr(1,k-1)+tpr(1,k)]*[fpr(1,k)-fpr(1,k-1)]/2;
end
overallauc=sum(area);
end
          
