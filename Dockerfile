FROM anaconda/miniconda:latest
RUN apt-get update && apt-get install -y ca-certificates && update-ca-certificates
# Limit conda channels to those that are free for use in any context
RUN echo -e "channel_priority: flexible\n\
channels:\n\
  - https://conda.anaconda.org/bioconda\n\
  - https://conda.anaconda.org/conda-forge\n\
  - https://conda.anaconda.org/r" \
  > /opt/miniconda3/.condarc
# install nextflow
RUN conda create -n nextflow26 -y python=3.12 nextflow=26
# preinstall all env that are otherwise built by nf upon running
COPY conda/*.yaml /tmp/
RUN for yml in /tmp/*.yaml;do conda create -n $(basename $yml)_env -f $yml python=3.12;done
# install libgomp needed by Kraken
RUN conda install -c conda-forge libgomp
ENV LD_LIBRARY_PATH=/opt/miniconda3/lib
# Kraken needs libgomp to be available under the name libgomp.so.1
RUN cp /opt/miniconda3/lib/libgomp.so.1.0.0 /opt/miniconda3/lib/libgomp.so.1
# Add the scripts and Nextflow workflow files to the container
COPY scripts/ /opt/karmii/scripts/
COPY *.nf     /opt/karmii/
# Update the Nextflow workflow files to use the preinstalled conda environments and fix tar extraction behavior to be compatible with the container environment
RUN for yml in /tmp/*.yaml;do sed -i "s/conda\/$(basename $yml)/\/opt\/miniconda3\/envs\/$(basename $yml)_env/g" /opt/karmii/*.nf;done
RUN sed -i "s/tar zxvf/tar --no-same-owner -zxvf/g" /opt/karmii/*.nf 
ENV JAVA_HOME=/opt/miniconda3/envs/nextflow26/bin/java
ENV PATH="/usr/local/bin:/opt/miniconda3/envs/nextflow26/bin:${PATH}"
