FROM ubuntu:12.04
RUN apt-get update
RUN apt-get install -y wget g++ gfortran python-numpy
RUN apt-get upgrade -y
RUN apt-get install -y python-scipy pdb2pqr make
RUN apt-get install -y python-pip && pip install --index-url https://pypi.python.org/simple/ --upgrade pip
RUN pip install peptidebuilder
ENV ATTRACTDIR=/attract/bin
ENV ATTRACTTOOLS=/attract/tools
ENV PYTHONPATH=/usr/share/pdb2pqr
RUN wget http://www.attract.ph.tum.de/services/ATTRACT/downloads/spydersilk-0.03.tar.gz
RUN pip install spydersilk-0.03.tar.gz
#RUN wget http://www.attract.ph.tum.de/services/ATTRACT/attract.tgz
#RUN tar xvzf attract.tgz
ADD . /attract
RUN cd attract/bin; make clean; make all -j4; exit 0
RUN ls ./attract
RUN cd ../attract/gui; python -c 'import spyder, attractmodel'

