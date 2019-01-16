function [S] = NBCLDA_GN2(A,B,C,D,E,F)
%%
% Function NBCLDA_GN2 is used to obtain the prediction results of
% LncRNA-Disease associations by applying the Naïve Bayesian Theory to GN2
% A: adjacency matrix for the miRNA-lncRNA associations
% A(i,j)=1 means lncRNA i is related to miRNA j
% B: adjacency matrix for the miRNA-disease associations
% B(i,j)=1 means disease i is related to miRNA j
% C: adjacency matrix for the lncRNA-disease associations
% C(i,j)=1 means lncRNA i is related to disease j
% D: adjacency matrix for the gene-lncRNA associations
% D(i,j)=1 means lncRNA i is related to gene j
% E: adjacency matrix for the gene-disease associations
% E(i,j)=1 means disease i is related to gene j
% F: adjacency matrix for the miRNA-gene associations
% F(i,j)=1 means miRNA i is related to gene j

%%
[n1,m1]=size(A);
[n2,m2]=size(B);
[n4,m4]=size(D);
[cn,cm]=size(C);
cc=sum(sum(C));
pld=(m1*m2-cc)/cc;
S1=ones(m1,m2);
S2=ones(m1,m2);
S3=ones(m1,m2);
N1=0;
N2=0;
N3=0;
for k=1:n1 
    Temp1=A(k,:)'*B(k,:);         %
    index=find(F(:,1)==k);
    if ~isempty(index)
       for i=1:length(index)
           Tempg=D(F(index(i),2),:)'*E(F(index(i),2),:);
           Temp=Temp1.*Tempg;
           Temp1=Temp1-Temp;
       end
    end
    x1=length(nonzeros(sum(Temp1,1)~=0)); 
    y1=length(nonzeros(sum(Temp1,2)~=0)); 
    Temp11=C.*Temp1;
    N1=sum(sum(Temp11));
    S11=Temp1.*(pld*(N1+1))/(x1*y1-N1+1);
    S11(S11==0)=1;
    S1=S11.*S1;    
end

for k=1:n4 
    Temp1=D(k,:)'*E(k,:);         %
    index=find(F(:,2)==k);
    if ~isempty(index)
       for i=1:length(index)
           Tempm=A(F(index(i),1),:)'*B(F(index(i),1),:);
           Temp=Temp1.*Tempm;
           Temp1=Temp1-Temp;
       end
    end
    Temp11=C.*Temp1;
    x2=length(nonzeros(sum(Temp1,1)~=0)); 
    y2=length(nonzeros(sum(Temp1,2)~=0)); 
  
    N2=sum(sum(Temp11));
    S22=Temp1.*(pld*(N2+1))/(x2*y2-N2+1);
    S22(S22==0)=1;
    S2=S22.*S2;   
end

for k=1:length(F)  
    Temp1=A(F(k,1),:)'*B(F(k,1),:);         
    Temp2=D(F(k,2),:)'*E(F(k,2),:);
    Temp=Temp1.*Temp2;
    x3=length(nonzeros(sum(Temp,1)~=0)); 
    y3=length(nonzeros(sum(Temp,2)~=0)); 
    Temp33=C.*Temp;
    N3=sum(sum(Temp33));
    S33=Temp.*(pld*(N3+1))/(x3*y3-N3+1);
    S33(S33==0)=1;
    S3=S33.*S3;   
end

S=S1.*S2.*S3;
S=log(S)/100;


 
