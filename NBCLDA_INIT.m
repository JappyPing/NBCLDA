clear
clc
datapath = [pwd,filesep,'data',filesep];
addpath(datapath);

%% load data to perform NBCLDA
load Disease_SemSim
DS = Disease_SemSim;
load Gene_Disease
GD = Gene_Disease;
load Gene_LncRNA
GL = Gene_LncRNA;
load LncRNA_Disease
LD = LncRNA_Disease;
load MiRNA_Disease
MD = MiRNA_Disease;
load MiRNA_Gene
MG = MiRNA_Gene;
load MiRNA_LncRNA
ML = MiRNA_LncRNA;
load MiRNAGene
mg = MiRNAGene;
load LncRNADisease
ld = LncRNADisease;
%% 
%Evaluation of NBCLDA_GN2_SD by using LOOCV
NBCLDA_LOOCV(ML,MD,GL,GD,mg,ld,LD,DS);
overallauc = positiontooverallauc(LD,ld);