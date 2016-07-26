# Docker image for VirMet

Docker containers are a lightweight alternative to virtual machines (VMs). We
prepared an image pre-built with all necessary dependencies that can be
downloaded and run with little effort. Docker runs on Linux, Windows and Mac
OS X. Within any of these hosts, it will run an isolated container inside which
VirMet is installed together with its dependencies. Make sure you have enough
disk space (hundreds of gigabytes) to save databases.

### Install docker

Follow the instructions reported [here](https://www.docker.com/products/docker)
to install docker on your host machine.

Check your installation with `docker version` and `docker run hello-world`.

### Pull docker image for VirMet

Simply type `docker pull ozagordi/virmet`. This will download approximately 1GB
of data: it contains a Linux Ubuntu distribution with all the required software
pre-installed.

### Start a docker container

The command to give is the following (backslashes can be omitted if the command
is given on a single line)

    [user@host ~ ]$ docker run -it --name my_virmet \
    -v your_path_to_some_data:/home/ubuntu/miseq_runs \
    ozagordi/virmet /bin/bash

Let's see this command in detail.

The above command will drop you in a bash shell (`/bin/bash`) in a new container
derived from the image `ozagordi/virmet`. This contained will be named
`my_virmet` (you can choose the name you prefer) in interactive mode (`-it`),
and it will mount a volume found on your host. Let's see an example of how
we can use this mount option.

Substitute `your_path_to_some_data` to a global path where you are storing
sequencing data. For example, if you have Miseq runs in `/data/MiSeqOutput`

    [user@host ~ ]$ ls -1 /data/MiSeqOutput
    160719_M01274_0159_000000000-AMA6L/
    160715_M02081_0162_000000000-AL5TJ/
    ...

then use the option `-v /data/MiSeqOutput:/home/ubuntu/miseq_runs`. What comes
before the `:` is the global path on your host machine to the Miseq runs; this
will be found inside the container at the path `/home/ubuntu/miseq_runs`. This
is a convenient way to make data on your host machine available to the
container, _i.e._ to VirMet.

The official reference for `docker run` is
[here](https://docs.docker.com/engine/reference/run/).

### Once in the container

The prompt is now different, now we are root inside the container
(identified with an ID string) and the prompt reads `root@18ef5268a7a4:/home/ubuntu#`.

First of all one has to activate the conda environment with the correct
installation of Python, pandas and Biopython. The command is

    root@18ef5268a7a4:/home/ubuntu# source activate test-virmet

We can check that the host directory is correctly mounted

    root@18ef5268a7a4:/home/ubuntu# ls /home/ubuntu/miseq_runs/
    160719_M01274_0159_000000000-AMA6L/
    160715_M02081_0162_000000000-AL5TJ/
    ...

Good! The data directories on our host are mounted inside the container.
Let's now make a directory to store our results

    root@18ef5268a7a4:/home/ubuntu# mkdir results
    root@18ef5268a7a4:/home/ubuntu# cd results

The first time we need to `fetch` (download) and `index` all databases, see
the relevant [instructions](preparation.md).

Once the databases are ready, we can finally run our first `wolfpack` command.

    root@18ef5268a7a4:/home/ubuntu/results# virmet wolfpack --run /home/ubuntu/miseq_runs/160719_M01274_0159_000000000-AMA6L &

and wait for the results to appear.

### Some suggestions

#### Mount multiple volumes

The `-v` option of `docker run` allows mounting of multiple volumes. For example
one could mount the directory with sequencing data, one to store results and
one for the databases by repeating `-v` multiple times

    [user@host ~ ]$ docker run -it --name my_virmet \
    -v your_path_to_some_data:/home/ubuntu/miseq_runs \
    -v your_path_to_db:/data/virmet_databases \
    -v your_path_to_results:/home/ubuntu/results \
    ozagordi/virmet /bin/bash

This is convenient if one wants to have available all files both within the
container and outside of it. In alternative, one can copy files from the
container to the host with `docker cp`, see
[reference](https://docs.docker.com/engine/reference/commandline/cp/). Users are
encouraged to get familiar with Docker and find the strategy that suits them best.

#### Leaving a container and attaching it again

A container can be left with `ctrl + p`, `ctrl + q` (release `ctrl` in between) and
attached later with `docker attach my_virmet` (or any other name). It might be
necessary to give an extra press to `return`.
If we give `ctrl + d` inside the container this is left and stopped, so one has
to start it again with `docker start my_virmet` before attaching it.
