FROM centos:centos8 as base

LABEL AUTHOR="Maximilian Miller <mmiller@bromberglab.org>"

RUN yum makecache && yum -y install\
 wget\
 perl\
 bzip2\
 gcc\
 ncurses-devel\
 zlib-devel\
 bzip2-devel\
 xz-devel\
 make

# variables
ENV SAMTOOLS_VERSION=1.9\
 HTSLIB_VERSION=1.9\
 TOOLS=/home/TOOLS/tools\
 TOOL_NAME=annovar\
 TOOL_VERSION=current\
 TARBALL_LOCATION=http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/\
 TARBALL=annovar.latest.tar.gz
ENV TARBALL_FOLDER=$TOOL_NAME\
 DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION\
 PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# annovar
RUN wget -q $TARBALL_LOCATION/$TARBALL &&\
 tar xf $TARBALL --wildcards *pl &&\
 rm -rf $TARBALL &&\
 cd $TARBALL_FOLDER &&\
 mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin &&\
 cp *pl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R &&\
 mkdir /databases &&\
 ln -s /databases/ $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb &&\
 cd ../ &&\
 rm -rf $TARBALL_FOLDER

# samtools, htslib
RUN cd / &&\
 wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 &&\
 tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 &&\
 cd /samtools-${SAMTOOLS_VERSION}/ && ./configure && make &&\
 mv /samtools-${SAMTOOLS_VERSION}/samtools /bin/ &&\
 cd htslib-${HTSLIB_VERSION}/ && ./configure && make &&\
 mv htsfile libhts.so* tabix bgzip /bin &&\
 rm -rf /samtools*

# cleanup
RUN yum clean all

# environment
RUN source /etc/profile.d/lang.sh
ENV PATH=$PATH:$TOOLS/$TOOL_NAME/current/bin

# setup bio-node
LABEL bio-node=v1.0

# set environment variables
WORKDIR /app

# set app ENTRYPOINT
ENTRYPOINT []

# set app CMD
CMD []
