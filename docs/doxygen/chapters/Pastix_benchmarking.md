## Automated PaStiX benchmarking

This section presents the work  put in place to autpmatically benchmark the PaStiX solver on the Plafrim platform.
It is composed of several modules that will be listed here.

### Initialization of the benchmark

**benchmark.yml**:

The script is launched through a gitlab runner with the tag
*plafrim*. It will use the `${PASTIX_DIR}/.gitlab/benchmark.yml` file
that sets the environment variables such as: the plafrim partition
(sirocco or bora), the amount of nodes, or the amount of cores
requested. Once the job is configured, it runs the `run.sh` script.

**run.sh**:

`benchmark.yml` calls the script `${PASTIX_DIR}/tools/bench/plafrim/run.sh` that sets
the constraints for the slurm job ( SLURM_CONSTRAINTS, PASTIX_BUILD_OPTIONS, MPI_OPTIONS, etc.) and
 launches the next script in the PaStiX guix environment: `slurm.sh`.

**slurm.sh**:

`${PASTIX_DIR}/tools/bench/plafrim/slurm.sh` launches the slurm job, and outputs the
progression of the benchmark in the job log. Each slurm job runs a job defined in the `pastix_guix.sh` script.

**pastix_guix.sh**:

`${PASTIX_DIR}/tools/bench/pastix_guix.sh` compiles PaStiX, launches the benchmarks thanks to JUBE
and pushes the results to an ElasticSearch database.

### JUBE

JUBE helps performing and analyzing benchmarks in a systematic way. It allows custom workflows to be
able to adapt to new architectures.

For each benchmark application the benchmark data is written out in a given format that enables JUBE
to deduct the desired information. This data can be parsed by automatic pre- and post-processing
scripts that draw information, and store it more densely for manual interpretation.

The JUBE benchmarking environment provides a script based framework to easily create benchmark sets,
run those sets on different computer systems and evaluate the results.

1- <u>Parameters</u>

First, we have to create an XML file that will contains the parameters for our runs.
It's structured like this:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<jube>
    <parameterset name="name_param">
        <parameter name="A"       type="string">A, B, C</parameter>
        <parameter name="N"       type="string">1, 2, 3</parameter>
        ...
        <parameter name="command" type="string">$A $N</parameter>
    </parameterset>
</jube>
```

The parameter "command" will be generated thanks to a Cartesian product:

```sh
> echo $command
A 1
A 2
A 3
B 1
B 2
B 3
C 1
C 2
C 3
```

2- <u>Run the benchmark</u>

Create a second XML file to launch the benchmark similar to this:

```xml
<benchmark name="benchmark" outpath="results">
    <step name="run_1" tag="run1">
        <use from="parameters.xml">name_param</use>
        <do>$command</do>
    </step>
    <step name="run_2" tag="run2">
        <use from="parameters.xml">name_param</use>
        <do>$command</do>
    </step>
    ...
</benchmark>
```

Then, launch the runs thanks to the following command. The tag is the field that allow us to recognize the run to do.

```sh
jube run run_file.xml --tag run1 run2 --include-path our/path/to/parameters.xml --id $id
```

The runs outputs are stored in the `results/id` directory. If `id` is not set in the command line,
a new one is automatically created.

3 - <u>Analyze the patterns</u>

JUBE can analyze the output files from the runs. To do that, a confifuration XML file tells JUBE how to extract the information from the output through regular expressions:

```xml
<jube>
    <patternset name="name_pattern">
        <pattern name="PATTERN1" type="string">regex</pattern>
        ...
    </patternset>
</jube>
```

Then, we can set the analyze command thanks to the following structure in our `run_file.xml`:

```xml
<analyser name="analyse">
    <!-- use a pattern set -->
    <use from="our/path/to/patterns.xml">name_pattern</use>
    <analyse step="run_1" tag="run1">
        <!-- file which should be scanned in results/id -->
        <file>stdout</file>
    </analyse>
    ...
</analyser>
```

And then call the analyze step thanks to:

```sh
jube analyse results
```

If the `id` is not set in the command line, the last one is used.

4- <u>Output the analyzed results</u>

Finally to extract the results in a CSV file from the run, the `run_file.xml` is used:

```xml
<result>
    <use>analyse</use> <!-- use existing analyser -->
    <table name="result" style="csv">
        <column>PATTERN1</column>
        ...
    </table>
</result>
```

And with the following command the result is converted:

```sh
jube result results --id $id > results.csv
```

## Send the results to the ElasticSearch database

The last step of the `pastix_guix.sh` script consists in sending the results to the ElasticSearch database. To perform this step, a python script parses the results from the runs to convert them to the ElasticSearch format, and then updates the database with the following command:

```sh
python3 ${PASTIX_DIR}/tools/bench/jube/add_result.py \
-e https://elasticsearch.bordeaux.inria.fr \
-t hiepacs -p "pastix" results.csv
```

The `-t` and `-p` options allow us to push the result as the index: *hiepacs-pastix_perf*

The rest of the work is done online, on the Kibana server.

## Visualize the runs on Kibana

### Index management

Once data has been posted with the *hiepacs-pastix_perf* index, we can visualize all the fields that have
been pushed to the ElasticSearch database.

Thanks to the index name, it is possible to create and a new index pattern  that will be used to create the tables and visualisations.

### Discovery table

In the discover page, we can select the wanted index pattern, and filter only the fields that
interest us as a **Search**. There are currently 2 **Search** used for the visualizations:
*pastix-perf_FR* and *pastix-perf_LR*.

### Visualization

Now we can visualize our data. First we decide to link it with one of our search and then we
will decide to the axis of the visualization.

For the moment, the X-axis corresponds to the pushed date and the Y-axis to the average values of the
selected data to analyze for this date.

### Dashboard

Here you can show a set of selected visualizations to have the wanted figures. Do to so, you have
to create/edit a Dashboard and you need to select the visualizations to add/remove/resize.

It seems that we can export a Dashboard or a single visualization as a PDF if we subscribe to the
ElasticSearch suite.
