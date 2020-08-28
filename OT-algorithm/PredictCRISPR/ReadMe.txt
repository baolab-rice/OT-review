This is the python version source codes for the CRISPR/Cas9 system off-target site prediction

*********************************************************************************************
1. About the files
   a. There is a source code file that contains the suffix of '.py':
      OfftargetPredict.py -- the main code for the prediction
      
   The code is programmed with python 2.7. Python packages such as numpy, scipy, sklearn, svmutil
   are required to be installed before running this file.
   
   b. The file folder "predict_results" saves the output results of the prediction in a csv file.
   c. The file folder "SVM_Model" contains the trained svm models (totally 40 models)
   d. The file minmaxs is a data file for our prediction model
   e. The file folder "example_input_files" gives two example input files where:
      offline_input.txt -- is the file that obtained from the offline software 'Cas-OFFinder'
      online_input.txt -- is an example file obtained from the website 'http://www.rgenome.net/cas-offinder/'

**********************************************************************************************
2. How to use this tool
   
    2.1 preparing the input files
       we use the genome wide candidate off-target site sequences of a given sgRNA on-target site sequence as
       the input of the tool. Users can prepare the input candidate off-target site sequences by two ways:
       A. from the web tool Cas-OFFinder: http://www.rgenome.net/cas-offinder/
          parameters: 

          a. ## PAM Type  ###############################
          The PAM Type should be " SpCas9 from Streptococcus pyogenes: 5'-NGG-3' " or 
          " SpCas9 from Streptococcus pyogenes: 5'-NRG-3' (R = A or G) ";
          
          b. ## Target Genome ############################
          Human and Mouse genomes can be selected (GRCh38/hg38 or hg19 for Human and mm10 for Mouse)
          
          c. ## Query Sequences ##########################
          a 20nt sequence of the on-target site should be entered into the blank;
          Mismatch Number is suggested to be set no bigger than 6 (<=6)
          DNA Bulge Size and RNA Bulge Size are remain to be 0
       
       With all the parameters being filled, then submit the job and wait and download the output file

       B. from the offline tool Cas-OFFinder which can be downloaded from: http://www.rgenome.net/cas-offinder/portable
       The detail of how to use this tool can be found from the website:  http://www.rgenome.net/cas-offinder/portable

   2.2 The Libsvm tool
       We applied the Libsvm version 3.22 to implement the binary classification.
       This tool can be downloaded from https://www.csie.ntu.edu.tw/~cjlin/libsvm/
       please install a python version of libsvm first before running our codes

   2.3 Run the tool
      The steps for running the tool are:
      a. Install the required packages such as numpy, scipy, sklearn, svmutil;
      b. preparing a candidate off-target site file according to previous 2.1;
      c. enter the folder containing the source code file "OfftargetPredict.py" under a command line tool;
      d. paste the following command:
         
         python OfftargetPredict.py [sgRNA] [input file type] [input file path]

      There are three required inputs:
         sgRNA --- string: is the 23nt on-target site sequence (20nt protospacer + 3nt PAM )
         type --- integer: to be 1 or 2, where 1 means the input candidate off-target sites are obtained with the online web tool;
                           2 means the candidate off-target sites are obtained with the offline software.
         filepath --- string: the complete filepath of the previous prepared candidate off-target site sequence
      
       e. find your output results under the folder xxx/predict_results/ (xxx means the current folder) 
  

      There is one output file:
         predict_results.csv --- the predicted scores that show the prabability of a candidate site to be a real off-target site
                                 0.5 is the threshold for label the candidate sites, where score>0.5 will be labeled as '1', otherwise '0',
                                 the predicted off-target sites have been ranked in the descending order of "score" field


####################################################################################################################################
####################################################################################################################################

Example command (under the working filefold that contain the OfftargetPredict.py):
    
       example1: predicting the off-target sites of the sgRNA with the spacer: 
                   sgSeq = 'AGGCACCGAGCTGTATGGTGTGG'
                   filepath = 'example_input_files\offline_input.txt'
                   type = 2
                   
                   python OfftargetPredict.py AGGCACCGAGCTGTATGGTGTGG 2 example_input_files\offline_input.txt
                   

       example2: predicting the off-target sites of the sgRNA with the spacer:
                   sgSeq = 'GACCCCCTCCACCCCGCCTCCGG'
                   filepath = 'example_input_files\online_input.txt'
                   type = 1
                   
                   python OfftargetPredict.py GACCCCCTCCACCCCGCCTCCGG 1 example_input_files\online_input.txt


##############################################################################################################################################
##############################################################################################################################################
Please contact Hui Peng: Hui.Peng-2@student.uts.edu.au if you encounter some problems when run the codes.

       
