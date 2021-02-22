#################################################################
# Dockerfile to build je and Trimmomatic
# images
# Based on Ubuntu
#  $ cd MintChIP
#  $ VERSION=dev
#  $ docker build -t jianhong/mintchip:$VERSION .  ## --no-cache
#  $ docker images jianhong/mintchip:$VERSION
#  $ docker push jianhong/mintchip:$VERSION
#  $ docker tag jianhong/mintchip:$VERSION jianhong/mintchip:latest
#  $ docker push jianhong/mintchip:latest
#  $ docker system prune -a
#  $ cd ~
#  $ docker pull jianhong/mintchip:latest
#  $ mkdir tmp4mintchip
#  $ docker run -it --rm -v ${PWD}/tmp4mintchip:/work jianhong/mintchip:latest
##################################################################
FROM ubuntu:20.10
LABEL authors="Jianhong Ou" \
      description="Docker image containing all software requirements for the jianhong/MintChIP pipeline"

ENV DEBIAN_FRONTEND="noninteractive" TZ="America/New_York"

# Install software from source
RUN cd ~ && \
    apt-get update --fix-missing && \
    apt-get install --yes wget python3 unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN ln -s python3 /usr/bin/python

RUN cd ~ && wget https://raw.githubusercontent.com/gbcs-embl/Je/master/dist/je_2.0.RC.tar.gz && \
    tar -xf je_2.0.RC.tar.gz && cd je_2.0.RC && \
    sed -i "s/bin\/sh/usr\/bin\/env bash/" je && \
    cp * /usr/local/sbin/ && cd .. && rm -rf je*

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && cd Trimmomatic-0.39 && \
    wget https://raw.githubusercontent.com/jianhong/MintChIP/dev/modules/trimmomatic/trimmomatic && \
    chmod +x trimmomatic && cp -r * /usr/local/sbin/ && \
    ln -s trimmomatic-0.39.jar /usr/local/sbin/trimmomatic.jar && \
    cd .. && rm -rf Trimmomatic*

WORKDIR /work
