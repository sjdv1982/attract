FROM ubuntu:12.04
RUN apt-get update
RUN apt-get install -y wget g++ gfortran python-numpy
RUN apt-get upgrade -y
RUN apt-get install -y python-scipy pdb2pqr make python-pip
RUN pip install  --index-url=https://pypi.python.org/simple/ peptidebuilder
ADD . /attract
ADD ./spydersilk-0.03 /spydersilk-0.03
RUN tar czf spydersilk-0.03.tar.gz spydersilk-0.03 && pip install spydersilk-0.03.tar.gz
RUN cd attract/bin; make clean; make all -j4; exit 0
RUN ls ./attract
RUN ln -s /usr/lib/gcc/x86_64-linux-gnu/4.6/libgfortran.so  /usr/lib/gcc/x86_64-linux-gnu/4.6/libgfortran.so.4
ENV LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/4.6
ENV ATTRACTDIR=/attract/bin
ENV ATTRACTTOOLS=/attract/tools
ENV PYTHONPATH=/usr/share/pdb2pqr
RUN cd attract/gui; python -c 'import spyder, attractmodel'
