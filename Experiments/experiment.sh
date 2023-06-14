#!/bin/bash -l
#SBATCH --job-name=SSTAT
#SBATCH --output=SSTAT.out
#SBATCH --error=SSTAT.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G



#1) ARGS="Nitta2015 Homo_sapiens";; #10
#2) ARGS="Jolma2013 Homo_sapiens";; #708
#3) ARGS="Jolma2013 Mus_musculus";; #134
#4) ARGS="Jolma2015 Homo_sapiens";; #658
#5) ARGS="Yin2017 Homo_sapiens";; #1794
#6) ARGS="fromYimeng pfm_composite_new";; #503
#7) ARGS="fromYimeng pfm_spacing_new";; #202




#*********** 3. sstat

#Returns the similarity between PFMs
# most parameters are the same as for cstat.

TF_PATH=../../TFBS/PWMs/Jolma2013/transfac/Homo_sapiens/
ls -d $TF_PATH* > ../data/all.list

TF_PATH=../../TFBS/PWMs/Jolma2013/transfac/Mus_musculus/
ls -d $TF_PATH* >> ../data/all.list

TF_PATH=../../TFBS/PWMs/Jolma2015/transfac/Homo_sapiens/
ls -d $TF_PATH* >> ../data/all.list

TF_PATH=../../TFBS/PWMs/Nitta2015/transfac/Homo_sapiens/
ls -d $TF_PATH* >> ../data/all.list #this took 7s and 3.5 MB

TF_PATH=../../TFBS/PWMs/Yin2017/transfac/Homo_sapiens/
ls -d $TF_PATH* >> ../data/all.list

TF_PATH=../../TFBS/PWMs/fromYimeng/transfac/pfm_composite_new/
ls -d $TF_PATH* >> ../data/all.list

TF_PATH=../../TFBS/PWMs/fromYimeng/transfac/pfm_spacing_new/
ls -d $TF_PATH* >> ../data/all.list

#Jolma2013: 50% GC-content, pseudocount regularisation, type I threshold 0.01, these give the
#same overall density of binding site occurrences for both PWMs tested and limit the effects
#of stringency and low affinity sites on the similarity score

#Yin2017 stringent type I error threshold 0.01 to limit the effect of low affinity similarities
#50 % GC-content background model
#pseudocount regularisation. This approach generally gives similar results as other common methods but performs better when Two
#otherwise dissimilar motifs share a common part

srun ../sstat .5 list:../data/all.list typeI 0.01 > ../results/all.txt

seff $SLURM_JOBID


#list:data/matrix.list
#./sstat <gc> <[list:]transfac-file> <threshold-method> <threshold-parameter> [<partial-execution>] [<return diagonal>]

#<gc>: gc content, e.g. '.4', for the background model
#<[file:]transfac-file>: file describing the position count matrices (PCM) in transfac format, e.g. 'data/A1.mat'.

#<threshold-method>: the treshold method:
  #typeI: set threshold such that typeI error is equal to threshold-parameter.
  #typeII: set threshold such that typeII error is equal to threshold-parameter.
  #balanced: set threshold such that typeI error equals typeII, threshold-parameter can be any number (is not used but has to be passed as parameter)
  #typeIext: set threshold to balanced threshold if it's possible such that typeI error is less or equal than threshold-parameter.
  #Otherwise, set typeI error equal to threshold-parameter.
  #threshold: threshold-parameter contains the threshold.
#[regularize]: if not set or set to a true value (1), the regularization method from Rahmann et al. 2005 is applied.
#Otherwise, we just add pseudocounts.
#[parameter|lambda|theta|rate|cpd]: if this parameter exist, the running time of the statistical calculation
#(without reading input/preparing PSSM, and so on) is printed. Depending on the choice of this parameter, following output is printed:
  #parameter: only xi, xi', xi'0, alpha, theta1
  #lambda: in addition: lambda1, lambda2
  #theta: in addition, all thetas until < precision
  #rate: in addition, rate r (give sequence length in next parameter) without theta
  #cpd: in addition, all P(X>=x) until < precision (give sequence length)
#[sequence-length]: length of the sequence

#The typeI error is measured as the probability of at least one
#false positive in a region of length 500 (see Pape et al. 2006: alpha_500).


#New parameters:

#<partial-execution>: integer i: if not given, whole similarity matrix is computed.
#if given, only the ith and the n-i th line of the similarity matrix are
#computed and return in special format (to be read by scluster). if -1 then simstat uses SGE cluster itself.

#<return diagonal>: default: 0 (false). If set to 1, we also return
#the similarity of each matrix and itself.
#(Useful for computing the variance for the univariate count distributions.)

#Output:

#Matrix with following columns:

#matrixA: first matrix
#matrixB: second matrix
#Smax: Similarities summarized by using the maximum
#Ssum: Similarities summarized by using the sum
#imax: Position with the highest similarity (B is shifted against fixed A!)
#bimaxp: maximum similarity is a reverse complementary hit (1) otherwise (0)
#alphaA: probability of a false positive for matrixA
#alphaB: probability of a false positive for matrixB

#Example:

#./sstat .4 list:data/matrix.list balanced .1

#Compute similarities between all pairs of matrices from data/matrix.list using a balanced threshold.



#*********** 8. Parallel Computing

#Some of the programs support parallel computing. Since we support OMP (one memory, multiple processors) and the Sun Grid Engine (multiple memory, multiple processors), we divide this section correspondingly.

#8a. OMP

#If you have an OMP ready compiler, you just have to uncomment the two lines in the Makefile:

#compileoption += -fopenmp
#linkoption = -lgomp -o

#And, perhaps, change the parameter for using openmp fitting to the compiler you use. We are using G++ compiler v. 4.2.0 for 64bit machines.

#If you have compiled the program with OMP enabled, the clustering will perform much faster in recomputing the similarities of each new representatives with all other nodes.

#8b. The Sun Grid Engine

#As the parameters suggested (above), the programs support the Sun Grid Engine - although we have to admit that the implementation is rather proprietary. Anyways, some inspection of the sge.h and sge.cpp should clarify the implementation and give the possibility to extend it. In fact, all classes which can use the SGE engine (CSimilarityMatrix and CClusterMatrix since it is inherited from CSimilarityMatrix but does not need any further adjustment.) are derived from the interface ISGEClient. Two adjustment might be needed:

#a) In the client class CSimilarityMatrix you might want to modify the path of sstat. We assume that it is contained in the path - then - you don't need any modifications.

#b) In sge.cpp, three main task are done - and might need some adjustments:

#(i) Initialization: The constructor of CSGEMaster needs a temporary directory (default: sgetemp). Be aware that each construction might delete files within this directory.

#(ii) Job Submission: Implemented in the method submit. Here, you have to change the format/directives and the program to submit the job (default:submit2sge) as well as the queue and other parameters such that they fit you environment.

#(iii) Waiting: After the jobs are submitted, the class waits until all jobs are done (method finish()). Here, we use the program qstat to see which jobs are done (using the job id caught at submission) and also perform some basic error handling. Depending on your output of qstat, your error logs and so on, you might want to change some code there, as well. (By the way, if we see that a job was finished successfully, we call the callback function sge_merge (for which we in fact use the interface ISGEClient) to read the output).

#*********** 9. Website & Comments

#Check out the website http://mosta.molgen.mpg.de for new revisions
#and an online version of the program. Send comments, suggestions etc. to utz.pape@molgen.mpg.de.

#*********** 10. Citation

#if you use the count statistic in any published work, please cite the following article (and, if possible, give a link to http://mosta.molgen.mpg.de):

#@UNPUBLISHED{pape2007a,
  #author = {Utz J. Pape and Sven Rahmann and Fengzhu Sun and Martin Vingron},
  #title = {{Compound Poisson approximation of DNA motif counts on both strands}},
  #note = {submitted to JCB},
  #year = {2007}
#}

#if you use the similarity or the clustering in any published work, please cite the following article (and, if possible, give a link to http://mosta.molgen.mpg.de):

#@UNPUBLISHED{pape2007b,
#  author = {Utz J. Pape and Sven Rahmann and Martin Vingron},
#  title = {{Natural Similarity Measures between Position Frequency Matrices with an Application to Clustering}},
#  note = {submitted to Bioinformatics},
#  year = {2007}
#}
