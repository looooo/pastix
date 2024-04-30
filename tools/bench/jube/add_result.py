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
import click
import csv
import json
import math
import sys
import time
from elasticsearch import Elasticsearch
from git import Repo
from typing import Any, Dict, List, Union

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


def format_entry(row: Row, commit_pastix: Repo, guix_commits: Dict) -> Dict[str, Any]:
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

def format_entry_stats(row: Row, commit_pastix: Repo, commit_stats: str,
                       guix_commits: Dict, es: Elasticsearch, es_index: str):
    """"format a result and compute stats: mean and stdev of gflop"""
    err = 0

    # format measures entry
    result = format_entry(row, commit_pastix, guix_commits)

    # prepare default gflops stats entry if not existing
    result_stats = {
        "Commit_date_pastix": result['Commit_date_pastix'],
        "Commit_sha_pastix": result['Commit_sha_pastix'],
        "Commit_sha_guix": result['Commit_sha_guix'],
        "Commit_sha_guix_hpc": result['Commit_sha_guix_hpc'],
        "Commit_sha_guix_hpcnonfree": result['Commit_sha_guix_hpcnonfree'],
        "Hostname": result['Hostname'],
        "Algorithm": result['Algorithm'],
        "Scheduler": result['Scheduler'],
        "NbMPI": result['NbMPI'],
        "NbGPU": result['NbGPU'],
        "CompressWhen": result['CompressWhen'],
        "SplitStart": result['SplitStart'],
        "mean": format(result['Gflops'], '.1f'),
        "stdev": format(result['Gflops']*0.1, '.1f')
    }

    if commit_stats != 'null':
        # search stats data for this commit and given input parameters
        search_param = {
        "query": {
          "bool": {
            "must": [
                { "match": { "Commit_sha_pastix": commit_stats }},
                { "match": { "Hostname" : result['Hostname'] }},
                { "match": { "Algorithm" : result['Algorithm'] }},
                { "match": { "Scheduler" : result['Scheduler'] }},
                { "match": { "NbMPI" : result['NbMPI'] }},
                { "match": { "NbGPU" : result['NbGPU'] }},
                { "match": { "CompressWhen" : result['CompressWhen'] }},
                { "match": { "SplitStart" : result['SplitStart'] }}
            ]
          }
        },
        "size": 1,
        "_source": ["mean", "stdev"],
        }
        response = es.search(index=es_index, body=search_param)
        elastic_docs2 = response["hits"]["hits"]

        if len(elastic_docs2) > 0:
            last_stats_data = elastic_docs2[0]["_source"]
            #print("last_stats_data ", last_stats_data)

            # compute formula from https://public.kitware.com/Wiki/CDash:Design#Test_Timing
            alpha = 0.3
            multiplier = 3

            previousMean = float(last_stats_data['mean'])
            previousSD = float(last_stats_data['stdev'])
            #print("previousMean ", previousMean)

            currentV = result['Gflops']
            # just to test: apply a perturbation
            #pert = random.uniform(-previousSD, previousSD)
            #currentV = result['Gflops'] + pert
            #print("currentV ", currentV)

            newMean = (1-alpha)*previousMean + alpha*currentV
            newSD = math.sqrt((1-alpha)*previousSD*previousSD + alpha*(currentV-newMean)*(currentV-newMean))

            # prepare stats data to put in database newMean and newSD
            result_stats = {
                "Commit_date_pastix": result['Commit_date_pastix'],
                "Commit_sha_pastix": result['Commit_sha_pastix'],
                "Commit_sha_guix": result['Commit_sha_guix'],
                "Commit_sha_guix_hpc": result['Commit_sha_guix_hpc'],
                "Commit_sha_guix_hpcnonfree": result['Commit_sha_guix_hpcnonfree'],
                "Hostname": result['Hostname'],
                "Algorithm": result['Algorithm'],
                "Scheduler": result['Scheduler'],
                "NbMPI": result['NbMPI'],
                "NbGPU": result['NbGPU'],
                "CompressWhen": result['CompressWhen'],
                "SplitStart": result['SplitStart'],
                "mean": format(newMean, '.1f'),
                "stdev": format(newSD, '.1f')
            }

            # check for regression
            thresholdSD = 0.1*previousMean
            if previousSD < thresholdSD:
                previousSD = thresholdSD
            maxAcceptableDiff = multiplier*previousSD
            diff = abs(currentV-previousMean)
            if diff > maxAcceptableDiff:
                print("Regression: inputs %(Hostname)s, %(Algorithm)s, %(Scheduler)s, %(NbMPI)s, %(NbGPU)s, %(CompressWhen)s, %(SplitStart)s " % result)
                print("Regression: outputs Gflops={0}, previousMean={1}, diff={2}, maxAcceptableDiff={3}".format(currentV, previousMean, diff, maxAcceptableDiff))
                err = 1

    return [result_stats, err]

@click.command()
@click.option("-d", "--directory", default=".", help="git working directory")
@click.option("-e", "--elastic-url", default="http://localhost:9200", help="elasticsearch instance url")
@click.option("-t", "--team", required=True, help="team name")
@click.option("-p", "--project", required=True, help="project name")
@click.argument("csv-files", nargs=-1)
def main(
    directory: str,
    elastic_url: str,
    team: str,
    project: str,
    csv_files: str,
):
    """Add a result to an elasticsearch database."""
    es = Elasticsearch(elastic_url)
    es_index = team + "-" + project + "_" + "perf"
    if not es.indices.exists(es_index):
        es.indices.create(es_index)

    # call this if mapping must be changed (e.g. add a new field)
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

    # load data from csv file
    requests = [
        request
        for file in csv_files
            for request in map(
                lambda row: format_entry(row, commit_pastix, guix_commits),
                open_csv(file)
            )
    ]

    # insert measures in database
    for request in requests:
        es.index(index=es_index.lower(), body=request)


    # compute stats: mean and stdev of gflops measured
    # database for stats
    es_index_stats = team + "-" + project + "_" + "stats"
    if not es.indices.exists(es_index_stats):
        es.indices.create(es_index_stats)

    # call this if mapping must be changed (e.g. add a new field)
    mapping_input_stats = {
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
            "FactoTime" : {"type" : "float"},
            "SolveTime" : {"type" : "float"},
            "RefineTime" : {"type" : "float"},
            "OrderTime" : {"type" : "float"},
            "SymbolTime" : {"type" : "float"},
            "ReorderTime" : {"type" : "float"},
            "BlendTime" : {"type" : "float"},
            "AnalyzeTime" : {"type" : "float"},
            "mean": {"type" : "float"},
            "stdev": {"type" : "float"}
        }
    }
    es.indices.put_mapping(index=es_index_stats, body=mapping_input_stats)

    # search last commit of database stats
    search_param = {
      "sort": [{"Commit_date_pastix": {"order": "desc"}}],
      "size": 1,
      "_source": ["Commit_sha_pastix"],
    }
    response = es.search(index=es_index_stats, body=search_param)
    elastic_docs = response["hits"]["hits"]
    last_stats_commit = 'null'
    if len(elastic_docs) > 0:
        # search last records commit
        last_stats_commit = elastic_docs[0]["_source"]['Commit_sha_pastix']
        print("Regression: mean and stdev taken from commit ", last_stats_commit)

    for file in csv_files:
        csvfile = open_csv(file)
        err = 0
        for row in csvfile:
            #print("row ", row)
            [entry, err2] = format_entry_stats(row, commit_pastix, last_stats_commit,
                                               guix_commits, es, es_index_stats)
            err = max(err, err2)
            # insert updated mean and stdev in database
            #print("entry ", entry)
            es.index(index=es_index_stats.lower(), body=entry)
        if err == 1:
            sys.exit(1)


if __name__ == "__main__":
    main()
