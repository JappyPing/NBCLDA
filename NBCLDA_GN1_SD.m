function [final_score]=NBCLDA_GN1_SD(A,B,C,D)
%%
% A: adjacency matrix for the lncRNA-miRNA associations
% A(i,j)=1 means lncRNA i is related to miRNA j
% B: adjacency matrix for the disease-miRNA associations
% B(i,j)=1 means disease i is related to miRNA j
% C: adjacency matrix for the lncRNA-disease associations
% C(i,j)=1 means lncRNA i is related to disease j
% D represents the disease semantic similarity
[S] = NBCLDA_GN1(A,B,C);
final_score = D*S';
