# Off-target scores calculation:  
The processing scripts adapted from previous publications: 

**SCRIPT: Function_based_algorithms.py**

Python 2 and pandas are required to run the python script.

The SCRIPT contains five function-based algorithms from previous published studies:

+ CCTOP off-target scores:
Stemmer, Manuel, et al. "CCTop: an intuitive, flexible and reliable CRISPR/Cas9 target prediction tool." PloS one 10.4 (2015).

+ Hsu score: CRISPOR review, Haeussler, Maximilian, et al. "Evaluation of off-target and on-target scoring algorithms and integration into the guide RNA selection tool CRISPOR." Genome biology 17.1 (2016): 148.  
https://github.com/maximilianh/crisporPaper

+ MIT(website) score: CRISPOR review, Haeussler, Maximilian, et al. "Evaluation of off-target and on-target scoring algorithms and integration into the guide RNA selection tool CRISPOR." Genome biology 17.1 (2016): 148.  
https://github.com/maximilianh/crisporPaper

+ CROP-IT score: CRISPOR review, Haeussler, Maximilian, et al. "Evaluation of off-target and on-target scoring algorithms and integration into the guide RNA selection tool CRISPOR." Genome biology 17.1 (2016): 148.  
https://github.com/maximilianh/crisporPaper

+ COSMID score: Cradick, Thomas J., et al. "COSMID: a web-based tool for identifying and validating CRISPR/Cas off-target sites." Molecular Therapy-Nucleic Acids 3 (2014): e214.  
https://crispr.bme.gatech.edu/

**Folder: CFDScoring**

Code for the CFD score, obtained from the authors.  

Doench, John G., et al. "Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9." Nature biotechnology 34.2 (2016): 184.  
Code also available in supplementary: https://www.nature.com/articles/nbt.3437#MOESM9 

**Folder: CNN_std**

Linked to the original github page of CNN_std algorithm. Can be installed based on authors instructions.

Lin, Jiecong, and Ka-Chun Wong. "Off-target predictions in CRISPR-Cas9 gene editing using deep learning." Bioinformatics 34.17 (2018): i656-i663.  
https://github.com/MichaelLinn/off_target_prediction

**Folder: CRISTA**

Code for the CRISTA score, obtained from the authors.  

Abadi, Shiran, et al. "A machine learning approach for predicting CRISPR-Cas9 cleavage efficiencies and patterns underlying its mechanism of action." PLoS computational biology 13.10 (2017): e1005807.  
http://crista.tau.ac.il/download.html

**Folder: Elevation**

Linked to the original github page of Elevation algorithm. Can be installed based on authors instructions. 

Listgarten, Jennifer, et al. "Prediction of off-target activities for the end-to-end design of CRISPR guide RNAs." Nature biomedical engineering 2.1 (2018): 38-47.  
https://github.com/microsoft/Elevation 

**Folder: PredictCRISPR**

Linked to the original github page of Elevation algorithm. Can be installed based on authors instructions. 

Peng, Hui, et al. "Recognition of CRISPR/Cas9 off-target sites through ensemble learning of uneven mismatch distributions." Bioinformatics 34.17 (2018): i757-i765.  
https://github.com/penn-hui/OfftargetPredict 



Default models were used for all machine learning algorithms.