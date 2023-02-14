FROM rocker/shiny-verse

RUN apt update && apt-get install -y --no-install-recommends libglpk-dev liblzma-dev libbz2-dev tcl-dev tk-dev libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

COPY  cran_pkg.R .
RUN Rscript cran_pkg.R

COPY  bioc_pkg.R .
RUN Rscript bioc_pkg.R

COPY  git_pkg.R .
RUN Rscript git_pkg.R

COPY ./scgenes/. /srv/shiny-server/.
RUN sudo chown -R shiny:shiny /srv/shiny-server/
