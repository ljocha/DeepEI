FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND=noninteractive 
ENV TZ=Europe/Prague

RUN apt update --fix-missing
RUN apt install -y python3-pip

RUN apt install -y python3-rdkit

RUN pip3 install scipy
RUN pip3 install pandas
RUN pip3 install sklearn
RUN pip3 install tqdm
RUN pip3 install jpype1

RUN apt install -y git
RUN pip3 install git+https://github.com/hcji/pycdk@master

RUN apt install -y default-jdk-headless

RUN apt install  -y r-base


#RUN Rscript -e 'install.packages("rcdklibs")'
#RUN Rscript -e 'install.packages("rJava")'
RUN Rscript -e 'install.packages("rcdk")'

RUN pip3 install rpy2
RUN pip3 install matchms

ENV JAVA_HOME=/usr/lib/jvm/default-java

RUN bash -c 'echo $JAVA_HOME/lib/server >/etc/ld.so.conf.d/java-for-R.conf'
RUN bash -c 'echo $JAVA_HOME/lib >>/etc/ld.so.conf.d/java-for-R.conf'
RUN ldconfig
