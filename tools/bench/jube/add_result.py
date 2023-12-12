"""
 @file add_result.py

 Python script to push benchamrk results to ElasticSearch database.

 @copyright 2021-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.3.2
 @author Tony Delarue
 @author Florent Pruvost
 @date 2023-07-21
"""
#!/usr/bin/env python3

from typing import Any, Dict, List, Union
from copy import deepcopy
import json
import click
import csv
import time
from git import Repo
from elasticsearch import Elasticsearch

Row = Dict[str, Union[str, float]]

Sched_str = [ "Sequential", "Static", "PaRSEC", "StarPU", "Dynamic" ]

def open_csv(filename: str) -> List[Dict[str, str]]:
    """
    Open a csv file a return it as dictionary.
    First row is titles.
    """
    csv_rows = []
    with open(filename) as csv_data:
        reader = csv.DictReader(csv_data)
        titles = reader.fieldnames
        for row in reader:
            csv_rows.append(
                {
                    title: row[title]
                    for title in titles
                }
            )
    return csv_rows


def format_entry(row: Row, mpivendor : str, commit_pastix: Repo, guix_commits: Dict) -> Dict[str, Any]:
    """"format a result"""
    commit_date_pastix = time.strftime( "%Y-%m-%d", time.gmtime(time.time()) ) + " 00:00:00"
    commit_sha_pastix  = str(commit_pastix.hexsha)
    hostname  = str( row.pop('hostname')  )
    algorithm = str( row.pop('algorithm') )
    scheduler = Sched_str[ int( row.pop('scheduler') ) ]
    nmpi  = int( row.pop('nmpi') )
    ngpu  = int( row.pop('ngpu') )
    when  = str( row.pop('when') )
    split = str( row.pop('SPLIT') )

    gflops = float( row.pop('GFLOPS_avg') )
    ftime  = float( row.pop('FTIME_avg')  )
    stime  = float( row.pop('STIME_avg')  )
    rtime  = float( row.pop('RTIME_avg')  )

    otime  = float( row.pop('OTIME')  )
    sytime = float( row.pop('SYTIME') )
    rotime = float( row.pop('ROTIME') )
    btime  = float( row.pop('BTIME')  )
    atime  = float( row.pop('ATIME')  )

    result = {
        "Commit_date_pastix": commit_date_pastix,
        "Commit_sha_pastix": commit_sha_pastix,
        "Commit_sha_guix": guix_commits["guix"],
        "Commit_sha_guix_hpc": guix_commits["guix-hpc"],
        "Commit_sha_guix_hpcnonfree": guix_commits["guix-hpc-non-free"],
        "Hostname": hostname,
        "Algorithm": algorithm,
        "Scheduler": scheduler,
        "NbMPI": nmpi,
        "NbGPU": ngpu,
        "CompressWhen": when,
        "SplitStart": split,
        "Gflops": gflops,
        "FactoTime": ftime,
        "SolveTime": stime,
        "RefineTime": rtime,
        "OrderTime": otime,
        "SymbolTime": sytime,
        "ReorderTime": rotime,
        "BlendTime" : btime,
        "AnalyzeTime": atime
    }
    return result


@click.command()
@click.option("-d", "--directory", default=".", help="git working directory")
@click.option("-e", "--elastic-url", default="http://localhost:9200", help="elasticsearch instance url")
@click.option("-t", "--team", required=True, help="team name")
@click.option("-p", "--project", required=True, help="project name")
@click.option("-h", "--host", required=True, help="host name")
@click.option("-m", "--mpi", required=True, help="MPI vendor (openmpi, nmad)")
@click.argument("csv-files", nargs=-1)
def main(
    directory: str,
    elastic_url: str,
    team: str,
    project: str,
    host: str,
    mpi: str,
    csv_files: str,
):
    """Add a result to an elasticsearch database."""
    es = Elasticsearch(elastic_url)
    es_index = team + "-" + project + "_" + "perf"
    if not es.indices.exists(es_index):
        es.indices.create(es_index)

    mapping_input = {
        "properties": {
            "Commit_date_pastix": {"type": "date", "format": "yyyy-MM-dd' 'HH:mm:ss"},
            "Commit_sha_pastix": {"type": "keyword"},
            "Commit_sha_guix": {"type": "keyword"},
            "Commit_sha_guix_hpc": {"type": "keyword"},
            "Commit_sha_guix_hpcnonfree": {"type": "keyword"},
            "Hostname": {"type": "keyword"},
            "Algorithm": {"type": "keyword"},
            "Scheduler": {"type": "keyword"},
            "NbMPI" : {"type" : "integer"},
            "NbGPU" : {"type" : "integer"},
            "CompressWhen" : {"type" : "keyword"},
            "SplitStart" : {"type" : "keyword"},
            "Gflops": {"type" : "float"},
            "FactoTime" : {"type" : "float"},
            "SolveTime" : {"type" : "float"},
            "RefineTime" : {"type" : "float"},
            "OrderTime" : {"type" : "float"},
            "SymbolTime" : {"type" : "float"},
            "ReorderTime" : {"type" : "float"},
            "BlendTime" : {"type" : "float"},
            "AnalyzeTime" : {"type" : "float"}
        }
    }

    es.indices.put_mapping(index=es_index, body=mapping_input)

    repo = Repo(directory, search_parent_directories=True)
    commit_pastix = repo.head.commit

    guix_commits = {"guix" : "",
                    "guix-hpc": "",
                    "guix-hpc-non-free" : ""}
    # collect guix commits info
    with open('guix.json') as f:
        guix_describe = json.load(f)
    for index_guix in guix_describe:
        if index_guix["name"] in guix_commits.keys() :
            guix_commits[ index_guix["name"] ] = index_guix["commit"]

    requests = [
        request
        for file in csv_files
            for request in map(
                lambda row: format_entry(row, mpi, commit_pastix, guix_commits),
                open_csv(file)
            )
    ]
    for request in requests:
        es.index(index=es_index.lower(), body=request)

if __name__ == "__main__":
    main()
