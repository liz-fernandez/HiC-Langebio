---
layout: page
title: Theoretical and Practical HiC Workshop
subtitle: Prerequisites
minutes: 5
---

## Software Prerequisites

For the course you will need to install R and R studio (If your version is a bit old please update it):

To download R:

[https://www.r-project.org/](https://www.r-project.org/)

You need to follow the CRAN download link and choose your favorite server and your operating system (Linux, Mac or Windows). 

To download R Studio:

[https://www.rstudio.com/](https://www.rstudio.com/)

You will also need to install Docker: 

[https://docs.docker.com/install/](https://docs.docker.com/install/)

Choose your system on the left hand side menu and follow the installation instructions. 

Once installed you will need to install the following Docker image via the command line:

~~~ {.bash}
$ docker pull lizfernandez/hic-langebio
~~~
To check it is running correctly use the following command:

~~~ {.bash}
$ docker run -i -t lizfernandez/hic-langebio /bin/bash
~~~
Once the prompt opens type:

~~~ {.bash}
$ HiCUP-master/hicup
~~~
You should see the help for hicup. You can exit the docker image just typing:

~~~ {.bash}
$ exit
~~~
