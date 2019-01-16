function [final_score] = NBCLDA_GN2_SD(A,B,C,D,E,F,DS)
%%
% DS represents the disease semantic similarity
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
[S] = NBCLDA_GN2(A,B,C,D,E,F);
final_score =DS*S';