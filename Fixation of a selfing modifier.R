#Fixation probability and the average time to fixation of a selfing modifier allele
#Author: Kuangyi Xu

#parameters
N=100; #population size
r0=0.8; #background selfing rate already in the population
c=0; #pollen discounting
d=0; #inbreeding depression
r=0.1; #effect of modifier allele
h=0.5; #dominance

#selfing rate of the three genotypes
rAA=r0+r;
rAa=r0+h*r;
raa=r0;

p<-c(); #gene frequency of allele A
p_AA<-c(); #genotype frequency of AA;
p_Aa<-c(); #genotype frequency of Aa;
p_aa<-c(); #genotype frequency of aa;

n_lost=0; # numebr of lost events
n_fix=0; #numebr of fixation events
t_lost=0; #average time to lost
t_fix=0; #average time to fix
nrep=10^4; #number of replications

for (rep in 1:nrep) {
  #initial condition
  t=0;
  p_aa[1]=(N-1)/N;
  p_Aa[1]=1/N;
  p_AA[1]=0;
  
  while (TRUE){
    t=t+1;
    if (p_aa[t]==1){
      nlost=nlost+1;
      t_lost=t_lost+t-1;
      break
    }
    if (p_AA[t]==1){
      nfix=nfix+1;
      t_fix=t_fix+t-1;
      break
    }
    
    p[t]=p_AA[t]+p_Aa[t]/2; 
    
    #outcrossing part
    #pollen pool
    pAm=p_AA[t]*(1-c*rAA) + p_Aa[t]*(1-c*rAa)/2;
    pam=p_Aa[t]*(1-c*rAa)/2 + p_aa[t]*(1-c*raa);
    pm=pAm+pam;
    p_pollen=pAm/pm;
    
    p_AA_out=p_AA[t]*(1-rAA)*p_pollen + p_Aa[t]*(1-rAa)*p_pollen/2;
    p_Aa_out=p_AA[t]*(1-rAA)*(1-p_pollen) + p_Aa[t]*(1-rAa)*((1-p_pollen)/2 + p_pollen/2)+ p_aa[t]*(1-raa)*p_pollen;
    p_aa_out=p_Aa[t]*(1-rAa)*(1-p_pollen)/2 + p_aa[t]*(1-raa)*(1-p_pollen);
    
    #selfing part
    p_AA_self=p_AA[t]*rAA + p_Aa[t]*rAa/4;
    p_Aa_self=p_Aa[t]*rAa/2;
    p_aa_self=p_Aa[t]*rAa/4 + p_aa[t]*raa;
    
    #proportion of selfed offspring
    p_self=p_AA_self+p_Aa_self+p_aa_self;
    
    #expected genotype frequency after reproduction and selection;
    W=(1-d)*p_self+(1-p_self);
    p_AAt=(p_AA_self*(1-d)+p_AA_out)/W;
    p_Aat=(p_Aa_self*(1-d)+p_Aa_out)/W;
    p_aat=(p_aa_self*(1-d)+p_aa_out)/W;
    
    #Sample the next generation
    draw=rmultinom(1, N, prob =c(p_AAt,p_Aat,p_aat))/N;

    #update genotype frequency
    p_AA[t+1]=draw[1];
    p_Aa[t+1]=draw[2];
    p_aa[t+1]=draw[3];
  }
}

p_fix=nfix/nrep; #fixation probability
t_fix=t_fix/nfix;#average time to fixation 
t_lost=t_lost/nlost; #average time to lost
