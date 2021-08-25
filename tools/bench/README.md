# PaStiX benchmarks

The scripts located here are used to analyze PaStiX performances in
a systematic way.  The workflow is as follows:

* Gitlab-ci
  We use a *schedule* job that will be triggered with a chosen
  frequency, the schedule job frequency can be tuned in the CI/CD part
  of the web interface
* Guix
  [Guix](https://guix.gnu.org/) is responsible for building an isolated and reproducible
  environment to build and execute PaStiX, Jube and some python
  scripts. We can use it on PlaFRIM because Guix and [Guix-HPC](https://gitlab.inria.fr/guix-hpc/guix-hpc-non-free)
  are installed.
* Jube
  [Jube](https://apps.fz-juelich.de/jsc/jube/jube2/docu/index.html) is used to drive the
  execution with different parameter spaces and to parse the results in csv files.
  Apart from the execution parameters such as the problem sizes, the number of ressources used
  and so on we also save in the database the commit date of pastix
  and the commit ids of pastix and guix channels to properly identify the software versions.
* Elasticsearch
  [Elasticsearch](https://www.elastic.co/fr/) is the database framework. The server is
  https://elasticsearch.bordeaux.inria.fr. It is only accessible from
  Inria's networks for now.
* Kibana
  [Kibana](https://www.elastic.co/fr/) is a web server to visualize the performances on graphs. It
  looks for data imported in the elasticsearch database. We want to be
  able to analyze the performances for each commit for which the
  scheduled job has been performed and to monitor some performances in
  the course of time/commits. Kibana server is hosted [here](https://kibana.bordeaux.inria.fr).
  It is only accessible from Inria's networks for now.
