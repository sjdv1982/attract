FROM ubuntu:18.04
RUN apt-get update && apt-get install -y nvidia-cuda-toolkit libboost-coroutine-dev libboost-system-dev libboost-program-options-dev
RUN apt-get install -y wget g++ gfortran python-numpy python-scipy make python-pip bc
RUN pip install  --index-url=https://pypi.python.org/simple/ --upgrade pip
RUN wget --no-check-certificate https://files.pythonhosted.org/packages/ff/f4/0ce39bebcbb0ff619426f2bbe86e60bc549ace318c5a9113ae480ab2adc7/biopython-1.76.tar.gz  && pip install biopython-1.76.tar.gz
RUN pip install  --index-url=https://pypi.python.org/simple/ peptidebuilder
RUN apt-get install -y pdb2pqr --fix-missing
ADD . /attract
ADD ./spydersilk-0.03 /spydersilk-0.03
RUN tar czf spydersilk-0.03.tar.gz spydersilk-0.03 && pip install spydersilk-0.03.tar.gz
ENV CUDADIR=usr/lib/nvidia-cuda-toolkit
RUN cd attract/bin; make cleanall; make all -j8
RUN ls /attract/bin/attract && chmod -R a+r+w+x /attract
ENV LD_LIBRARY_PATH=/attract/bin:/usr/lib/nvidia-cuda-toolkit/lib64:/usr/lib/nvidia-cuda-toolkit/lib64
ENV ATTRACTDIR=/attract/bin
ENV ATTRACTTOOLS=/attract/tools
ENV PYTHONPATH=/usr/share/pdb2pqr
RUN cd attract/gui; python2 -c 'import spyder, attractmodel'
