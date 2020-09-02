# Bioinformatics analysis for the Off-Target (OT) review

## Data collection:

**Folder: data_original**  
The original data sources of our manually curated true positive list collected from 9 studies (as described in Supplementary Table S3): 

+ 2019_Gomez_Ospina  
Gomez-Ospina, Natalia, et al. "Human genome-edited hematopoietic stem cells phenotypically correct Mucopolysaccharidosis type I." Nature communications 10.1 (2019): 1-14. 

+ 2020_Vaidyanathan  
Vaidyanathan, Sriram, et al. "High-efficiency, selection-free gene repair in airway stem cells from cystic fibrosis patients rescues cftr function in differentiated epithelia." Cell stem cell 26.2 (2020): 161-171. 

+ 2013_Fu  
Fu, Yanfang, et al. "High-frequency off-target mutagenesis induced by CRISPR-Cas nucleases in human cells." Nature biotechnology 31.9 (2013): 822-826. 

+ 2014_Cho  
Cho, Seung Woo, et al. "Analysis of off-target effects of CRISPR/Cas-derived RNA-guided endonucleases and nickases." Genome research 24.1 (2014): 132-141. 

+ 2019_Pavel-Dinu  
Pavel-Dinu, Mara, et al. "Gene correction for SCID-X1 in long-term hematopoietic stem cells." Nature communications 10.1 (2019): 1-15.

+ 2019_Park  
Park, So Hyun, et al. "Highly efficient editing of the beta-globin gene in patient-derived hematopoietic stem and progenitor cells to treat sickle cell disease." Nucleic Acids Research 47.15 (2019): 7955-7972. 

+ 2015_Kim  
Kim, Daesik, et al. "Digenome-seq: genome-wide profiling of CRISPR-Cas9 off-target effects in human cells." Nature methods 12.3 (2015): 237.

+ 2015_Wang  
Wang, Xiaoling, et al. "Unbiased detection of off-target cleavage by CRISPR-Cas9 and TALENs using integrase-defective lentiviral vectors." Nature biotechnology 33.2 (2015): 175. 

+ 2016_Kim  
Kim, Daesik, et al. "Genome-wide target specificities of CRISPR-Cas9 nucleases revealed by multiplex Digenome-seq." Genome research 26.3 (2016): 406-415. 



For each gRNA, off-target sites were screened by Cas-Offinder(http://www.rgenome.net/cas-offinder/) allowing up to 4 mismatches and 1 base DNA/RNA bulge. 


## Off-target scores calculation:  
**Folder: OT-algorithm**  
The processing scripts adapted from previous publications (as described in Table 3). Please see details of each algorithm inside the folder: 

+ CCTOP off-target scores were computed based on the formula in the original paper:  
Stemmer, Manuel, et al. "CCTop: an intuitive, flexible and reliable CRISPR/Cas9 target prediction tool." PloS one 10.4 (2015).

+ Code for the MIT score (Hsu score), and CROP-IT score was adapted from the CRISPOR review:  
Haeussler, Maximilian, et al. "Evaluation of off-target and on-target scoring algorithms and integration into the guide RNA selection tool CRISPOR." Genome biology 17.1 (2016): 148.  
https://github.com/maximilianh/crisporPaper

+ Code for the CFD score was obtained from the authors:  
Doench, John G., et al. "Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9." Nature biotechnology 34.2 (2016): 184.  
Code available in supplementary: https://www.nature.com/articles/nbt.3437#MOESM9 

+ elevation algorithm was installed based on authors instructions:  
Listgarten, Jennifer, et al. "Prediction of off-target activities for the end-to-end design of CRISPR guide RNAs." Nature biomedical engineering 2.1 (2018): 38-47.  
https://github.com/microsoft/Elevation 

+ predictCRISPR algorithm was installed based on authors instructions:  
Peng, Hui, et al. "Recognition of CRISPR/Cas9 off-target sites through ensemble learning of uneven mismatch distributions." Bioinformatics 34.17 (2018): i757-i765.  
https://github.com/penn-hui/OfftargetPredict 

+ CNN_std algorithm was installed based on authors instructions:  
Lin, Jiecong, and Ka-Chun Wong. "Off-target predictions in CRISPR-Cas9 gene editing using deep learning." Bioinformatics 34.17 (2018): i656-i663.  
https://github.com/MichaelLinn/off_target_prediction

+ Code for the CRISTA score was obtained from the authors:  
Abadi, Shiran, et al. "A machine learning approach for predicting CRISPR-Cas9 cleavage efficiencies and patterns underlying its mechanism of action." PLoS computational biology 13.10 (2017): e1005807.  
http://crista.tau.ac.il/download.html

+ Code for the COSMID score was obtained from the authors:  
Cradick, Thomas J., et al. "COSMID: a web-based tool for identifying and validating CRISPR/Cas off-target sites." Molecular Therapy-Nucleic Acids 3 (2014): e214.  
https://crispr.bme.gatech.edu/

Default models were used for all machine learning algorithms.

## Data processing and visualization:
**Folder: data_processed**

+ Supplementary Table S2: The full off-target dataset used in the performance assessment. Raw data for generating Figure 2A and Figure S1. 

+ Supplementary Table S3: The off-target dataset of novel gRNAs used in the performance assessment.  Raw data for generating Figure 2B and Figure S2. Potential off-target sites and scores of 4 gRNAs that were not included in any training sets of machine-learning-based algorithms.

+ ROC_PRC.py: The script of making Figure S1 and S2.
Python 3 required. 
Dependencies can be installed using conda:  
  conda install numpy scikit-learn pandas scipy  
  conda install -c conda-forge matplotlib

