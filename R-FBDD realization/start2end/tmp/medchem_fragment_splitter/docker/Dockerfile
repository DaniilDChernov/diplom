# Ubuntu 22.04 contains fresh RDKit in its repos
FROM ubuntu:22.04

# Install only what is required for the project above the clean system
RUN apt-get update && apt-get upgrade -y && \
	apt-get install -y python3-rdkit  

# CharegeParameterFitting project source is mapped here in a container
VOLUME /opt/mcfs

# By default, build the project and run functional tests in it
CMD ["cd /opt/mcfs/files; python3 /opt/mcfs/files/test.py", "/opt/mcfs"]
