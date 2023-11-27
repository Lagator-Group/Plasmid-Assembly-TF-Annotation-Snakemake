#DeepTFactor
##Procedure

**Note**: 
This source code was developed in Linux, and has been tested in Ubuntu 16.04 with Python 3.6.

1. Clone the repository

        git clone https://bitbucket.org/kaistsystemsbiology/deeptfactor.git

2. Create and activate virtual environment

        conda env create -f environment.yml
        conda activate deeptfactor

3. To use GPU for the computation, install an appropriate version of pytorch (and cuda).  
    - Please refer to <https://pytorch.org>


##Example


- Run DeepTFactor

        python tf_running.py -i ./Dataset/example_tf.fasta -o ./result -g cpu
        python tf_running.py -i ./Dataset/example_tf.fasta -o ./result -g cuda:1