FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for IKMB Virus pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/virus-pipe-1.0/bin:/opt/biobloom/bin:/opt/spades/3.15.0/bin:$PATH

RUN apt-get update && apt-get -y install procps make gcc  git build-essential autotools-dev automake libsparsehash-dev libboost-all-dev \
cmake zlib1g-dev coreutils grep gawk

RUN cd /opt && \
        git clone https://github.com/simongog/sdsl-lite.git && \
        cd sdsl-lite && \
        ./install.sh /usr/local/

RUN cd /opt && \
	wget https://github.com/bcgsc/biobloom/releases/download/2.3.1/biobloomtools-2.3.1.tar.gz && \
	tar -xvf biobloomtools-2.3.1.tar.gz && rm biobloomtools-2.3.1.tar.gz && cd  biobloomtools-2.3.1 && \
        ./configure --prefix=/opt/biobloom && make install && \
        cd /opt && rm -Rf biobloomtools-2.3.1

RUN cd /opt && \
	mkdir -p spades && cd spades && \
	wget https://github.com/ablab/spades/releases/download/v3.15.0/SPAdes-3.15.0-Linux.tar.gz && \
	tar -xvf SPAdes-3.15.0-Linux.tar.gz && rm *.tar.gz && \
	mv SPAdes-3.15.0-Linux 3.15.0
	
