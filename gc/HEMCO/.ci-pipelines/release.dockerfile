FROM liambindle/penelope:0.1.0-ubuntu16.04-gcc7-netcdf4.5.0-netcdff4.4.4

# Make a directory to install HEMCO to
RUN mkdir /opt/hemco && mkdir /opt/hemco/bin

# Copy the HEMCO repository (".") to /hemco-src
# This means this docker build command's context must be 
# HEMCO's root source code directory
COPY . /hemco-src
RUN cd /hemco-src \
&&  mkdir build

# Commands to properly set up the environment inside the container
RUN echo "module load gcc/7" >> /init.rc \
&&  echo "spack load hdf5" >> /init.rc \
&&  echo "spack load netcdf" >> /init.rc \
&&  echo "spack load netcdf-fortran" >> /init.rc \
&&  echo "export PATH=$PATH:/opt/hemco/bin" >> /init.rc

# Make bash the default shell
SHELL ["/bin/bash", "-c"]

# Build Standard and copy the executable to /opt/hemco/bin
RUN cd /hemco-src/build \
&&  cmake -DRUNDIR=IGNORE -DCMAKE_COLOR_MAKEFILE=FALSE .. \
&&  make -j install \
&&  cp hemco_standalone /opt/hemco/bin/hemco-standard \
&& rm -rf /hemco-src/build/*

RUN rm -rf /hemco-src

RUN echo "#!/usr/bin/env bash" > /usr/bin/start-container.sh \
&&  echo ". /init.rc" >> /usr/bin/start-container.sh \
&&  echo 'if [ $# -gt 0 ]; then exec "$@"; else /bin/bash ; fi' >> /usr/bin/start-container.sh \
&&  chmod +x /usr/bin/start-container.sh
ENTRYPOINT ["start-container.sh"]
