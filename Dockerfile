# this is our first build stage, it will not persist in the final image
FROM python:3.11-slim-buster
WORKDIR /project

# install git
RUN apt-get update && \
    apt install -y gcc git r-base openssh-client wget && \
    rm -rf /var/lib/apt/lists/* && \
    R -e "install.packages(c('pacman'), dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    R -e "pacman::p_load(optparse, cli, conumee, minfi, parallel, tibble, tidyr, dplyr, data.table, gtools)" && \
    R -e "pacman::p_load_gh('crukci-bioinformatics/rascal')" 

# add credentials on build
ARG SSH_PRIVATE_KEY
RUN mkdir /root/.ssh/
RUN echo "${SSH_PRIVATE_KEY}" > /root/.ssh/id_rsa && \
    chmod 600 /root/.ssh/id_rsa && \
    touch /root/.ssh/known_hosts

RUN ssh-keyscan github.com >> /root/.ssh/known_hosts
RUN git clone git@github.com:pbousquets/beast-mcmc-flipflop.git

# Set up environment variables
ENV JAVA_HOME=/project/jdk-21.0.1
ENV PATH=$JAVA_HOME/bin:$PATH
ENV PATH=/project/beast-mcmc-flipflop/release/Linux/Phyfumv1.0_RC1/bin:$PATH

# install the correct version of java!!!!!
RUN wget -O jdk.tar.gz https://download.oracle.com/java/21/latest/jdk-21_linux-x64_bin.tar.gz  && tar zxvf jdk.tar.gz

# install Ant
RUN mkdir -p /opt/ant/ && \
    wget http://archive.apache.org/dist/ant/binaries/apache-ant-1.9.8-bin.tar.gz -P /opt/ant && \
    tar -xvzf /opt/ant/apache-ant-1.9.8-bin.tar.gz -C /opt/ant/ && \
    rm -f /opt/ant/apache-ant-1.9.8-bin.tar.gz

#Add Sonarqube lib
RUN wget http://downloads.sonarsource.com/plugins/org/codehaus/sonar-plugins/sonar-ant-task/2.3/sonar-ant-task-2.3.jar -P /opt/ant/apache-ant-1.9.8/lib/
#Set Ant Home in Docker container for Jenkins
ENV ANT_HOME=/opt/ant/apache-ant-1.9.8
#Set Ant Params in Docker for Jenkins
ENV ANT_OPTS="-Xms256M -Xmx512M"
#Update Path
ENV PATH="${PATH}:${HOME}/bin:${ANT_HOME}/bin"

# Install beast
RUN cd /project/beast-mcmc-flipflop && ant linux

RUN pip install phyfum

ENTRYPOINT ["phyfum", "run"]
