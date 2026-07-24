          seed = -1
seqfile = iqtree.dummy.phy
treefile = iqtree.rooted.nwk
       outfile = mcmctree.out

         ndata = 2
       seqtype = 2  * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge =   * safe constraint on root age, used if no fossil for root.

       BDparas = 1 1 0.1 M   * birth, death, sampling
   rgene_gamma = 2 20 1   * gamma prior for overall rates for genes
  sigma2_gamma = 1 10 1    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 0: .1  .1  .1  .1 .1 .1  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

*** These parameters control the MCMC run
***  Note: Total number of MCMC iterations will be burnin + (sampfreq * nsample)

         print = 1
        burnin = 10000
      sampfreq = 10
       nsample = 10000


*** The following parameters only needed to run MCMCtree with exact likelihood (usedata = 1)
*** no need to change anything for approximate likelihood (usedata = 2)

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5    * alpha for gamma rates at sites
         ncatG = 4    * No. categories in discrete gamma

     cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

*** Note: Make your window wider (100 columns) before running the program.
