FROM ubuntu:14.04
MAINTAINER osvaldo.zagordi@gmail.com

RUN apt-get update -y && apt-get install -y \
    build-essential \
    bwa \
    git \
    libncurses5-dev \
    libwww-perl \
    r-cran-ggplot2 \
    seqtk \
    tabix \
    unzip \
    wget \
    zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# install prinseq
WORKDIR /tmp
RUN echo "installing prinseq"
RUN wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz \
    -O /tmp/prinseq-lite-0.20.4.tar.gz \
    && tar -xvf prinseq-lite-0.20.4.tar.gz \
    && install -v prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin/prinseq \
    && rm -rf prinseq-lite-0.20.4*

# install blast main tools
WORKDIR /tmp
RUN echo "installing blast"
RUN wget -nv -O - ftp://ftp.ncbi.nlm.nih.gov/blast//executables/blast+/2.3.0/ncbi-blast-2.3.0+-x64-linux.tar.gz | \
    tar -xz
WORKDIR /tmp/ncbi-blast-2.3.0+/bin
RUN install makeblastdb blastn /usr/local/bin/
WORKDIR /home/ubuntu
RUN rm -rf /tmp/ncbi-blast-2.3.0+

# install samtools
WORKDIR /tmp
RUN echo "installing samtools"
RUN wget -nv https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 \
    -O samtools-1.3.tar.bz2
RUN tar xvfj samtools-1.3.tar.bz2
RUN cd samtools-1.3 && make && make install
RUN cd ..
RUN rm -rf samtools-1.3*

# install edirect
WORKDIR /usr/local
RUN echo "installing edirect"
RUN perl -MNet::FTP -e \
  '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login; \
   $ftp->binary; $ftp->get("/entrez/entrezdirect/edirect.zip");';
RUN unzip -u -q edirect.zip
RUN rm edirect.zip
RUN ./edirect/setup.sh \
    && /bin/bash -c "source $HOME/.bash_profile"

# now the python environment with conda

# Force bash always
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
WORKDIR /home/ubuntu
# maybe move CONDA_ENV_PATH to /opt/miniconda
ENV CONDA_ENV_PATH /opt/miniconda
ENV MY_CONDA_PY3ENV "python3.5"

# TODO: remove miniconda.sh
RUN wget -nv https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
&& bash miniconda.sh -b -p $CONDA_ENV_PATH
ENV PATH $CONDA_ENV_PATH/bin:$PATH
RUN hash -r \
    && conda config --set always_yes yes --set changeps1 no \
    && conda update -q conda
# if psutil is installed as it asks to do, an error is thrown
# RUN conda install psutil
# Useful for debugging any issues with conda
#RUN conda info -a

# install deps from environment.yml and create the env
COPY environment.yml ./
RUN conda env create -q python=3.5
# source activate does not work, see
# http://stackoverflow.com/questions/37945759/condas-source-activate-virtualenv-does-not-work-within-dockerfile
# so we manually adjust the path
ENV PATH /usr/local/edirect/:/opt/miniconda/envs/test-virmet/bin:$PATH
# RUN python -m coverage --version
RUN pip install codecov

# without setting the locale test_common fails due to non ascii code in downloaded files
RUN locale-gen "en_US.UTF-8"
ENV LC_ALL="en_US.UTF-8"

WORKDIR /opt
RUN git clone --depth=50 --branch=master https://github.com/ozagordi/VirMet.git \
    && cd /opt/VirMet/ \
    && python setup.py install

WORKDIR /home/ubuntu
ENTRYPOINT "/bin/bash"
