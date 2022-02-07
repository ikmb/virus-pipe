FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for IKMB Virus pipeline"

COPY environment.yml /

RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/virus-pipe-1.3/bin:/opt/biobloom/bin:/opt/spades/3.15.0/bin:/opt/vt:$PATH

RUN apt-get -y update && apt-get -y install procps make gcc  git build-essential autotools-dev automake libsparsehash-dev libboost-all-dev \
cmake zlib1g-dev coreutils librest-client-perl librole-rest-client-perl libjson-parse-perl libjson-perl

RUN cd /opt && \
        git clone https://github.com/simongog/sdsl-lite.git && \
        cd sdsl-lite && \
        ./install.sh /usr/local/

RUN cd /opt && \
	wget https://github.com/bcgsc/biobloom/releases/download/2.3.1/biobloomtools-2.3.1.tar.gz && \
	tar -xvf biobloomtools-2.3.1.tar.gz && rm biobloomtools-2.3.1.tar.gz && cd  biobloomtools-2.3.1 && \
        ./configure --prefix=/opt/biobloom && make install && \
        cd /opt && rm -Rf biobloomtools-2.3.1
	
RUN /opt/conda/envs/virus-pipe-1.3/bin/snpEff download NC_045512.2

RUN cd /opt && git clone https://github.com/atks/vt.git && cd vt \
	&& git checkout 0.577 && make
