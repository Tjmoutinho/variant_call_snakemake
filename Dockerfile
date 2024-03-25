FROM condaforge/miniforge3:latest
ENV DEBIAN_FRONTEND=noninteractive

# Install general deps
RUN apt-get update && apt-get install -y \
    git \
    curl \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libvcflib-tools \
    libvcflib-dev \
    libtool \
    wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install conda env from yaml
COPY env.yml /usr/src/vcf/env.yml
RUN mamba env create -f /usr/src/vcf/env.yml

RUN mamba clean --all -y

# Copy the entrypoint script into the container and set it to executable
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Clone and compile BWA
RUN git clone https://github.com/lh3/bwa.git /bwa
WORKDIR /bwa
RUN make
ENV PATH="/bwa:$PATH"

# Install Samtools
RUN mkdir /opt/samtools
RUN curl -L https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 -o samtools-1.19.2.tar.bz2 && \
    tar -xf samtools-1.19.2.tar.bz2 && \
    cd samtools-1.19.2 && \
    ./configure --prefix=/opt/samtools/ && \
    make && \
    make install
ENV PATH="${PATH}:/opt/samtools/bin"

ENTRYPOINT ["/entrypoint.sh"]

CMD ["/bin/bash"]