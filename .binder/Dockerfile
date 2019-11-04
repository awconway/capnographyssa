FROM rocker/binder:3.6.1
LABEL maintainer='Aaron Conway'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}

RUN installGithub.r rstudio/gt

RUN wget https://github.com/awconway/capnographyssa/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-07-10/')); devtools::install_deps()"

RUN rm DESCRIPTION.1; exit 0
