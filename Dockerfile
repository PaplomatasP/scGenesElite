FROM rocker/shiny-verse

RUN apt update && apt-get install -y --no-install-recommends libglpk-dev liblzma-dev libbz2-dev tcl-dev tk-dev libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

COPY scgenes/Install_Packages/ /opt/scgenes/Install_Packages/
WORKDIR /opt/scgenes
RUN Rscript Install_Packages/install_all.R

COPY ./scgenes/. /srv/shiny-server/.
RUN sudo chown -R shiny:shiny /srv/shiny-server/
WORKDIR /srv/shiny-server
