FROM ubuntu:12.04
RUN apt-get update && apt-get install -y wget g++ gfortran python-numpy python-scipy
RUN apt-get upgrade -y
RUN apt-get install -y make python-pip --fix-missing
RUN pip install  --index-url=https://pypi.python.org/simple/ --upgrade pip
RUN wget --no-check-certificate https://files.pythonhosted.org/packages/ff/f4/0ce39bebcbb0ff619426f2bbe86e60bc549ace318c5a9113ae480ab2adc7/biopython-1.76.tar.gz  && pip install biopython-1.76.tar.gz
RUN pip install  --index-url=https://pypi.python.org/simple/ peptidebuilder
RUN apt-get -y remove python-numpy python-scipy && pip install  --index-url=https://pypi.python.org/simple/ peptidebuilder numpy scipy
RUN apt-get install -y pdb2pqr --fix-missing
ADD . /attract
ADD ./spydersilk-0.03 /spydersilk-0.03
RUN tar czf spydersilk-0.03.tar.gz spydersilk-0.03 && pip install spydersilk-0.03.tar.gz
RUN cd attract/bin; make clean; make all -j4; exit 0
RUN ls /attract/bin/attract && chmod -R a+r+w+x /attract
RUN ln -s /usr/lib/gcc/x86_64-linux-gnu/4.6/libgfortran.so  /usr/lib/gcc/x86_64-linux-gnu/4.6/libgfortran.so.4
ENV LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/4.6
ENV ATTRACTDIR=/attract/bin
ENV ATTRACTTOOLS=/attract/tools
ENV PYTHONPATH=/usr/share/pdb2pqr
RUN cd attract/gui; python2 -c 'import spyder, attractmodel'
