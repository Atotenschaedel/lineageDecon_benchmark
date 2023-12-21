# -*- coding: utf-8 -*-
"""
Created on Tue Aug 91 2023

@author: aschedl
"""
import yaml
import numpy as np
import pandas as pd
from pathlib import Path

# filter result dataframes from prediction pipeline per tool
def filter_dataframe(df: pd.DataFrame, threshold, tool_name):
    assert "timepoint" in df.columns, f"no timepoint column in dataframe for {tool_name}"

    # change all columns to float except timepoint
    for col in df.columns:
        if col != "timepoint":
            try: 
                df[col] = df[col].astype(float)
            except: pass

    # add missing columns if necessary
    if not "others" in df:
        df["others"] = np.NaN
    else: pass
    if not "sample_name" in df:
        df["sample_name"] = (tool_name + "-") + df["timepoint"].astype(str)
    else: pass
    if not "tool_name" in df:
        df["tool_name"] = tool_name
    else: pass

    # filter dataset
    # sum of all lineages under threshold
    for index, row in df.iterrows():
        other_sum = 0
        for col in df.select_dtypes(include=['float64']).columns:
            if isinstance(row[col], float) and row[col] <= threshold and not np.isnan(row[col]): 
                other_sum = other_sum + row[col]
                df.loc[index, col] = np.NaN
        df.loc[index, "others"] = other_sum

    # add to 100% if below, add difference to "others"
    for index, row in df.iterrows():
        sum_total = 0
        for col in df.select_dtypes(include=['float64']).columns:
            if isinstance(row[col], float) and not np.isnan(row[col]):
                sum_total = sum_total + row[col]
        if sum_total < 1.0:
            diff = 1.0 - sum_total
            df.loc[index, "others"] = row["others"] + diff

    df = df.round(2)
    df.dropna(axis=1, how='all', inplace=True)
    return df

# get either children or parent lineage from variant list, return extended list
def get_all_relatives(variant_list: list, children_data: dict, relation: str):
    assert relation == "children" or relation == "parent", "relation needs to be either 'children' or 'parent'"

    # get number of refernce lineages + their parent/children
    variant_list_adjusted = variant_list.copy()

    for var in variant_list:
        try:
            variant_list_adjusted.extend(children_data[var][relation])
        except KeyError:
            if var != "other": print(f"{var} not found in pango json file.")

    return set(variant_list_adjusted)


def find_number_RefLineages(config: dict, children_data: dict, w_dir: Path):
# find number of lineage in referece dataset for specific tools + their children (adj.)
    number_lineages = {}

    for tool in config["TOOLS"]:
        if config["TOOLS"][tool]:
            tool_name = tool.lower()
            variant_list = []

            if tool_name == "vaquero":
                ref_file = w_dir /  "bin" / "VaQuERo" / "resources"
                ref_file = ref_file / ("mutations_list_grouped_pango_codonPhased_" + config["VAQUERO"]["MARKER_DATE"] + "_" + config["VAQUERO"]["MARKER_LOCATION"] +".csv")
                df = pd.read_csv(ref_file, sep=',')

                for index, row in df.iterrows():
                    l = row["Variants"].split(";")
                    for v in l:
                        variant_list.append(v)

                number_lineages[tool_name] = len(set(variant_list))
                number_lineages[tool_name+"_adj"] = len(get_all_relatives(variant_list, children_data, relation="parent"))

            elif tool_name == "freyja":
                ref_file = w_dir / "bin" / "freyja" / "lineages.yml"
                with open(ref_file, "r") as f:
                    df = yaml.safe_load(f)

                variant_list = []
                for l in range(len(df)):
                    variant_list.append(df[l]["name"])

                number_lineages[tool_name] = len(set(variant_list))
                number_lineages[tool_name+"_adj"] = len(get_all_relatives(variant_list, children_data, relation="parent"))

            elif tool_name == "lollipop":
                ref_file = w_dir /  "reference" / "GISAID" / "phe-genomics" / "phe2cojac" / "variants_pangolin.yaml"
                with open(ref_file, "r") as f:
                    df = yaml.safe_load(f)

                number_lineages[tool_name] = len(df["variants_pangolin"])
                number_lineages[tool_name+"_adj"] = len(get_all_relatives(list(df["variants_pangolin"].values()), children_data, relation="parent"))

            elif tool_name == "vlq":
                variant_list = []
                ref_file = w_dir /  "reference" / "GISAID" / "VLQ_reference_set" / "lineages.txt"
                with open(ref_file, "r") as f:
                    for line in f:
                        variant_list.append(line.split()[0])
                number_lineages[tool_name] = len(set(variant_list))
                number_lineages[tool_name+"_adj"] = len(get_all_relatives(variant_list, children_data, relation="parent"))
            
    # always add reference set for simulation        
    variant_list = config["SIMULATION"]["VARIANTS"].values()
    variant_list = [item[1].strip("-cons") for item in variant_list]

    number_lineages["simulation"] = len(config["SIMULATION"]["VARIANTS"])
    number_lineages["simulation_adj"] = len(get_all_relatives(variant_list, children_data, relation="parent"))

    return number_lineages


# count TP = true positive / TN = true negative / FN = false negative /  FP = false positve
def TrueFalseCounts(sim_hits: list, tool_hits: list, tool_name: str, number_of_refLineages: dict, children_data: dict):
    # 'others' column not counted
    if "others" in sim_hits:
        sim_hits.remove("others")
    if "others" in tool_hits:
        tool_hits.remove("others")

    # strict mode
    TP = 0
    FP = 0
    FN = 0

    for lineage in tool_hits:
        if lineage in sim_hits:
            TP += 1
        else:
            FP += 1
    
    for lineage in sim_hits:
        if lineage not in tool_hits:
            FN += 1

    TN = number_of_refLineages[tool_name] - (len(sim_hits))

    # adjusted - children allowed, no penalty
    TP_adj_child = 0
    FP_adj_child = 0
    FN_adj_child = 0

    for lineage in tool_hits:
        try:
            if lineage in sim_hits or children_data[lineage]["parent"] in sim_hits:
                TP_adj_child += 1
            else:
                FP_adj_child += 1
        except KeyError:
            FP_adj_child += 1

    for lineage in sim_hits:
        hit = False
        if lineage in tool_hits:
            hit = True
            break
        else:
            try:
                for child in children_data[lineage]["children"]:
                    if child in tool_hits:
                        hit = True
                        break
            except KeyError:
                pass
        if hit: 
            pass
        else:
            FN_adj_child += 1

    TN_adj_child = number_of_refLineages[tool_name+"_adj"] - (len(sim_hits))

    return [TP, TN, FP, FN, TP_adj_child, TN_adj_child, FP_adj_child, FN_adj_child]


# calcuate MSE (mean squared error) of relative abundance differences
def delta_relAbundance(sim_data: pd.Series, tool_data: pd.Series, children_data: dict):
    # both dataframe should only have 1 row each
    assert len(sim_data)==1, "Simulation Dataframe is more then 1 row long, relative abundance data needs to be measured row-wise"

    # clean up
    sim_data.reset_index(drop=True, inplace=True)

    tool_data= tool_data.to_frame().T
    tool_data.dropna(axis=1, how='all', inplace=True)
    tool_data.drop(columns=["sample_name", "tool_name"], inplace=True)
    tool_data.reset_index(drop=True, inplace=True)

    assert len(tool_data)==1, "Tool Dataframe is more then 1 row long, relative abundance data needs to be measured row-wise"

    # calculate MSE (mean squared error) - sum((o-e)^2) / n

    # strict mode
    MSEsum = {}
    for lineage in tool_data.columns:
        if lineage in sim_data.columns:
            obs = tool_data[lineage][0]
            exp = sim_data[lineage][0]
        else:
            obs = tool_data[lineage][0]
            exp= 0
        
        MSEsum[lineage] = (obs - exp) ** 2
    
    for lineage in sim_data.columns:
        if lineage not in tool_data.columns:
            obs = 0 
            exp= sim_data[lineage][0]
            MSEsum[lineage] = (obs - exp) ** 2

    # adjusted - children allowed, no penalty
    MSEsum_adj = {}

    # first collect all abundances from lineage + children
    for lineage in tool_data.columns:
        MSEsum_adj[lineage] = []
        #check for perfect hits
        if lineage in sim_data.columns:
            obs = tool_data[lineage][0]
            MSEsum_adj[lineage].append(obs)
        else:
            # check if pangolin file has children/ parent data for lineage
            got_data = False
            try:
                children_data[lineage]["parent"]
                got_data = True
            except KeyError:
                got_data = False
            
            # if parent of lineage is in simulated dataset, add to parent abundance
            # remove child lineage from count
            if got_data:
                if children_data[lineage]["parent"] in sim_data.columns:
                    obs = tool_data[lineage][0]
                    try: 
                        MSEsum_adj[children_data[lineage]["parent"]].append(obs)
                    except KeyError:
                        MSEsum_adj[children_data[lineage]["parent"]] = [obs]
                    MSEsum_adj.pop(lineage)

                else:
                    # lineage or parent not part of simulted dataset - false postive
                    # false positive - expected value is 0
                    obs = tool_data[lineage][0]
                    exp = 0
                    MSEsum_adj[lineage] = (obs - exp) ** 2
            else:
                obs = tool_data[lineage][0]
                exp = 0
                MSEsum_adj[lineage] = (obs - exp) ** 2

    #sum over all abundances and compare with simulation
    for lineage in MSEsum_adj.keys():
        # lineage is in simulated dataset - true postive
        if lineage in sim_data.columns:
            obs = sum(MSEsum_adj[lineage])
            exp = sim_data[lineage][0]
            MSEsum_adj[lineage] = (obs - exp) ** 2

    # look for false negative hits
    for lineage in sim_data.columns:
        if lineage != "others":
            try:
                children = children_data[lineage]["children"]
                children.append(lineage)
            except KeyError:
                children = [lineage]

            # if no intersection lineage and none of its children is in tools list
            intersection = list(set(children) & set(tool_data.columns))
            if len(intersection) == 0:
                # observed value is 0
                #print(f"could not find {lineage} or any children in ", tool_data.columns)
                obs = 0 
                exp= sim_data[lineage][0]
                MSEsum_adj[lineage] = (obs - exp) ** 2
        else:
            if "others" not in tool_data.columns:
                obs = 0
                exp = sim_data[lineage][0]
                MSEsum_adj[lineage] = (obs - exp) ** 2

    if len(MSEsum) == 0:
        MSEsum = np.NaN
    else:
        MSEsum = sum(MSEsum.values()) / len(MSEsum)

    if len(MSEsum_adj) == 0:
        MSEsum_adj = np.NaN
    else:
        MSEsum_adj = sum(MSEsum_adj.values()) / len(MSEsum_adj)

    return [MSEsum, MSEsum_adj]


def get_simulationStats(config:dict, data: pd.DataFrame, experiment_name:str, w_dir: Path):

    data["coverage_avg"] = np.NaN
    data["coverage_sd"] = np.NaN
    data["uniformity_wg_per"] = np.NaN
    data["MAPQ_avg"] = np.NaN

    for index, row in data.iterrows():
        sample_name = row["sample"]

        # get coverage of sample
        i_file = w_dir / "experiments" / experiment_name / "results"/ "variantCall" / "00_stats" / sample_name/(sample_name + "_cov.tsv")
        df = pd.read_csv(i_file, sep="\t")

        # overall coverage
        data.loc[index, "coverage_avg"] = df[sample_name].mean()
        data.loc[index, "coverage_sd"] = df[sample_name].std()

        # uniformity overall - percentage of bases in all targeted regions (or whole‑genome) 
        # that is covered by at least X%.
        qc_passed = 0

        for i, r in df.iterrows():
            if int(r[sample_name]) >= int(config["POSTPRED"]["MIN_READ_COUNT"]):
                qc_passed += 1

        data.loc[index, "uniformity_wg_per"] = (qc_passed / len(df)) * 100

        # get average quality score of mapped reads
        i_file = w_dir / "experiments" / experiment_name / "results" / "variantCall"/ "00_stats" / sample_name / (sample_name + "_virus.stats")
        with open(i_file) as f:
            first_line =f.readline()
        data.loc[index, "MAPQ_avg"] = first_line.split("=")[1].strip()

    return data

