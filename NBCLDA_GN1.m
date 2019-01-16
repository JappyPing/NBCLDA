function [S1] = NBCLDA_GN1(A,B,C)
%%
%A: adjacency matrix for the lncRNA-miRNA associations
%A(i,j)=1 means lncRNA i is related to miRNA j
%B: adjacency matrix for the disease-miRNA associations
%B(i,j)=1 means disease i is related to miRNA j
%C: adjacency matrix for the lncRNA-disease associations
%C(i,j)=1 means lncRNA i is related to disease j
%%
[n1,m1] = size(A);
[n2,m2] = size(B);
cc = sum(sum(C));
pld = (m1*m2-cc)/cc;
S1 = ones(m1,m2);
N1 = 0;
S = ones(m1,m2);
for k = 1:n1 
    Temp = A(k,:)'*B(k,:);        
    Temp1 = C.*Temp;
    x1 = length(nonzeros(sum(Temp,1)~=0));
    y1 = length(nonzeros(sum(Temp,2)~=0));
    N1 = sum(sum(Temp1));
    S  = Temp.*(pld*(N1+1))/(x1*y1-N1+1);
    S(S==0)=1;
    S1 = S.*S1;
end
S1 = log(S1)/100;


           

