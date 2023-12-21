# -*- coding: utf-8 -*-
"""
Created on Fri July 29 2023

@author: aschedl
"""
# get parameters from config file
import yaml, os, json, argparse
import numpy as np
from pathlib import Path
import pandas as pd
import PostPred_functions as func

# surpress performance warning
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


def main(w_dir: Path, config: dict, exp_name: str, debug=False):
    # set experiment path
    exp_path = w_dir / 'experiment' / exp_name

    # create folders if needed
    if not os.path.exists((w_dir / "experiments" / exp_name / "results" / "postPrediction").absolute()):
        os.makedirs((w_dir / "experiments" / exp_name / "results" / "postPrediction").absolute())

    # set minimal lineage abundance threshold
    threshold = config['POSTPRED']['LINEAGE_MIN_THRESHOLD']

    # import simulated data + metadata
    ## read out metadata
    meta_df = pd.read_csv(w_dir / 'experiments' / exp_name / 'simulation' / (exp_name + '_metadata.tsv'), sep='\t')

    meta_df["timepoint"] = ""
    for index, row in meta_df.iterrows():
        meta_df.at[index, "timepoint"] = row["sample"].split("-")[1]

    meta_df["timepoint"] = meta_df["timepoint"].astype(int)
    meta_df["sample_date"] = pd.to_datetime(meta_df["sample_date"], format='%Y-%m-%d')

    ## read out sim data
    sim_df = pd.DataFrame({"timepoint": pd.Series(dtype=int)})

    for f in (w_dir / 'experiments' / exp_name / 'simulation' / 'abundances').iterdir():
        if "summary" not in f.stem and f.suffix == ".tsv":
            df = pd.read_csv(f, sep='\t', header=None, names=["lineage", "abundance"])

            # check if abundance is 0-100 or 0-1
            if df['abundance'].sum() == 100.0: 
                df['abundance'] = df['abundance'] / 100

            df['lineage'] = df['lineage'].str.replace('-cons','')
            df = df.T
            df.columns = df.iloc[0]
            df = df[1:]        
            df["timepoint"] = f.stem.split("-")[1]
            df["timepoint"] = df["timepoint"].astype(int)
            sim_df = pd.concat([sim_df, df], ignore_index=True).sort_values(by=["timepoint"])

    # change all columns to float except timepoint
    for col in sim_df.columns:
        if col != "timepoint":
            sim_df[col] = sim_df[col].astype(float)
            if not config["SIMULATION"]["USE_REAL_TIMECOURSE"]:
                sim_df[col] = sim_df[col] / 100
    

    sim_df["sample_name"] = "simul-" + sim_df["timepoint"].astype(str)
    sim_df["tool_name"] = "simulation"
    sim_df.replace(0.00, np.nan, inplace=True)
    sim_df = sim_df.merge(meta_df, on ="timepoint")

    if not debug:
        print()
        print("Simulation data: ")
        print(sim_df)

    # write summary file
    sim_df.to_csv(w_dir / "experiments" / exp_name / "results" / "postPrediction" / "simulation_summary.csv", index=False)
    merged_df = sim_df.copy()

    # read out freyja data
    if config["TOOLS"]["FREYJA"]:
        frey_df = pd.DataFrame({"timepoint": pd.Series(dtype=int)})

        if (w_dir / 'experiments' / exp_name / 'results' / 'freyja').is_dir():
            for f in (w_dir / 'experiments' / exp_name / 'results' / 'freyja').iterdir():
                if f.suffix == ".tsv":
                    names = []
                    values = []
                    with open(f) as file:
                        lines = file.readlines()
                        for line in lines:
                            line = line.replace('\n', "")
                            if line.startswith("lineages"):
                                line = line.replace("lineages\t", "")
                                names.extend(line.split(" "))
                            elif line.startswith("abundances"):
                                line = line.replace("abundances\t", "")
                                values.extend(line.split(" "))

                    df = pd.DataFrame([values], columns=names)
                    df["timepoint"] = f.stem.split("-")[1]
                    df["timepoint"] = df["timepoint"].astype(int)
                    frey_df = pd.concat([frey_df, df], ignore_index=True).sort_values(by=["timepoint"])

            # filter dataset
            frey_df = func.filter_dataframe(frey_df.copy(), threshold, tool_name='freyja')
        else:
            print("No Results found for freyja.")
            frey_df["timepoint"] = sim_df[["timepoint"]].copy()
            frey_df["tool_name"] = 'freyja'
            frey_df["sample_name"] = frey_df["tool_name"] + "-" + frey_df["timepoint"]
        
        # write to summary file
        frey_df.to_csv(w_dir / "experiments" / exp_name / "results" / "postPrediction" / "freyja_summary.csv", index=False)

        # merge with main df
        merged_df = pd.concat((merged_df, frey_df), axis = 0)

        if not debug:
            print()
            print("freyja data: ")
            print(frey_df)

            
    # read out lollipop
    if config["TOOLS"]["LOLLIPOP"]:
        if (w_dir / 'experiments' / exp_name / 'results' / 'lollipop' / 'deconvoluted.tsv').is_file():
            lol_df = pd.read_csv(w_dir / 'experiments' / exp_name / 'results' / 'lollipop' / 'deconvoluted.tsv', sep='\t')

            lol_df["date"] = pd.to_datetime(lol_df["date"], format='%Y-%m-%d')
            lol_df = lol_df.merge(meta_df, left_on ='date', right_on='sample_date')
            lol_df.drop(columns=["location", "date", "sample", "sample_date"], inplace=True)
            lol_df = lol_df.pivot(index='timepoint',columns='variant', values='proportion')
            lol_df.rename(columns = {"undetermined": "others"}, inplace=True)
            lol_df.reset_index(inplace=True)

            # filter dataset
            lol_df = func.filter_dataframe(lol_df.copy(), threshold, tool_name='lollipop')
        else:
            print("No Results found for lollipop.")
            lol_df["timepoint"] = sim_df[["timepoint"]].copy()
            lol_df["tool_name"] = 'lollipop'
            lol_df["sample_name"] = lol_df["tool_name"] + "-" + lol_df["timepoint"]

        # write to summary file
        lol_df.to_csv(w_dir / "experiments" / exp_name / "results" / "postPrediction" / "lollipop_summary.csv", index=False)

        # merge with main df
        merged_df = pd.concat((merged_df, lol_df), axis = 0)

        if not debug:
            print()
            print("lollipop data: ")
            print(lol_df)

    # read out VLQ
    if config["TOOLS"]["VLQ"]:
        VLQ_df = pd.DataFrame({'timepoint': pd.Series(dtype=int)})
        if (w_dir / 'experiments' / exp_name / 'results' / 'VLQ').is_dir:
            for f in (w_dir / 'experiments' / exp_name / 'results' / 'VLQ').iterdir():
                if f.is_file():
                    #colum names: [variant, tpm, freq(%), adj_freq(%)]
                    df = pd.read_csv(f, sep='\t', skiprows=3, header=None, names=["variant", "tpm", "freq(%)", "adj_freq(%)"])
                    df = df[df["adj_freq(%)"] != 0.0]
                    df.reset_index(inplace=True, drop=True)
                    df.drop(columns=["tpm", "freq(%)"], inplace=True)
                    
                    # adjust to abundance between 1-0
                    df["adj_freq(%)"] = df["adj_freq(%)"] / 100

                    df = df.T
                    df.columns = df.iloc[0]
                    df = df[1:]
                    df.reset_index(inplace=True, drop=True)
                    
                    df["timepoint"] = int(f.stem.split("-")[1].split("_")[0])
                    VLQ_df = pd.concat([VLQ_df, df], ignore_index=True).sort_values(by=["timepoint"])

            # filter dataset
            VLQ_df = func.filter_dataframe(VLQ_df.copy(), threshold , tool_name='VLQ')
        else:
            print("No Results found for VLQ.")
            VLQ_df["timepoint"] = sim_df[["timepoint"]].copy()
            VLQ_df["tool_name"] = 'VLQ'
            VLQ_df["sample_name"] = VLQ_df["tool_name"] + "-" + VLQ_df["timepoint"]
            
        # write to summary file
        VLQ_df.to_csv(w_dir / "experiments" / exp_name / "results" / "postPrediction" / "vlq_summary.csv", index=False)

        # merge with main df
        merged_df = pd.concat((merged_df, VLQ_df), axis = 0)

        if not debug:
            print()
            print("VLQ data: ")
            print(VLQ_df)

    # read out vaquero
    if config["TOOLS"]["VAQUERO"]:
        if (w_dir / 'experiments' / exp_name / 'results' / 'vaquero' / 'globalFittedData.csv').is_file():
            vaquero_df = pd.read_csv(w_dir / 'experiments' / exp_name / 'results' / 'vaquero' / 'globalFittedData.csv', sep='\t')
            vaquero_df["sample_date"] = pd.to_datetime(vaquero_df["sample_date"], format='%Y-%m-%d')
            vaquero_df = vaquero_df[vaquero_df["value"] != 0.0]
            vaquero_df = vaquero_df.merge(meta_df, on='sample_date')
            vaquero_df.drop(columns=["sample_id", "sample", "sample_date", "LocationID", "LocationName"], inplace=True)
            vaquero_df = vaquero_df.pivot(index='timepoint',columns='variant', values='value')
            vaquero_df.reset_index(inplace=True)

            # filter dataset
            vaquero_df = func.filter_dataframe(vaquero_df.copy(), threshold, tool_name='vaquero')
        else:
            print("No Results found for vaquero.")
            vaquero_df = sim_df[["timepoint"]].copy()
            vaquero_df["tool_name"] = 'vaquero'
            vaquero_df["sample_name"] = vaquero_df["tool_name"] + "-" + str(vaquero_df["timepoint"])

        # write to summary file
        vaquero_df.to_csv(w_dir / "experiments" / exp_name / "results" / "postPrediction" / "vaquero_summary.csv", index=False)

        # merge with main df
        merged_df = pd.concat((merged_df, vaquero_df), axis = 0)

        if not debug:
            print()
            print("vaquero data: ")
            print(vaquero_df)
            
    print()
    # merge all datasets
    merged_df.reset_index(inplace=True, drop=True)
    merged_df = merged_df.sort_values(by=["timepoint"])

    # data for adjusted counts - children allowed - data from cornelius römer repo
    with open(w_dir / "reference" / "pango-sequences" / "data" / "pango-consensus-sequences_summary.json") as json_file:
        json_data = json.load(json_file)

    number_of_refLineages = func.find_number_RefLineages(config, json_data, w_dir)
    if not debug:
        print("Number of lineages on reference set: ", number_of_refLineages)

    main_df = pd.DataFrame({'timepoint': pd.Series(dtype=int)})

    for ti_p in sim_df["timepoint"]:
        ti_p_df = merged_df.loc[merged_df['timepoint'] == ti_p].copy()
        ti_p_df.drop(columns=["sample", "sample_date", "timepoint"], inplace=True)
        ti_p_df.dropna(axis=1, how='all', inplace=True)
        ti_p_df.reset_index(inplace=True, drop=True)
        
        # copy simulated data for timepoint
        sim_ti_p = sim_df.loc[sim_df['timepoint'] == ti_p].copy()
        sim_ti_p.drop(columns=["sample", "sample_date", "sample_name", "tool_name", "timepoint"], inplace=True)
        sim_ti_p = sim_ti_p.loc[:, (sim_ti_p != 0).any(axis=0)]
        sim_ti_p.dropna(axis=1, how='all', inplace=True)
        
        ### TPR / FPR strict mode
        # TP = true positive
        # TN = true negative 
        # FN = false negative
        # FP = false positve

        ti_p_df["TP"] = np.NaN
        ti_p_df["TN"] = np.NaN
        ti_p_df["FP"] = np.NaN
        ti_p_df["FN"] = np.NaN

        ti_p_df["TP_adj_child"] = np.NaN
        ti_p_df["TN_adj_child"] = np.NaN
        ti_p_df["FP_adj_child"] = np.NaN
        ti_p_df["FN_adj_child"] = np.NaN

        ti_p_df["relAb_MSE"] = np.NaN
        ti_p_df["relAb_MSE_adj_child"] = np.NaN
        
        c = ti_p_df.columns.to_numpy()
        res = [c[x].tolist() for x in ti_p_df.notna().to_numpy()]
        
        for index, row in ti_p_df.iterrows():
            # only take columns not NA on this timepoint for this tool, remove tool_name, sample_name and others
            tool_data = res[index]
            tool_data.remove("sample_name")
            tool_data.remove("tool_name")

            # returns [TP, TN, FP, FN, TP_adj_child, TN_adj_child, FP_adj_child, FN_adj_child]
            count = func.TrueFalseCounts(sim_hits=list(sim_ti_p.columns), 
                                            tool_hits=tool_data, 
                                            tool_name=row["tool_name"].lower(), 
                                            number_of_refLineages=number_of_refLineages, 
                                            children_data=json_data)

            ti_p_df.loc[index, "TP"] = count[0]
            ti_p_df.loc[index, "TN"] = count[1]
            ti_p_df.loc[index, "FP"] = count[2]
            ti_p_df.loc[index, "FN"] = count[3]

            ti_p_df.loc[index, "TP_adj_child"] = count[4]
            ti_p_df.loc[index, "TN_adj_child"] = count[5]
            ti_p_df.loc[index, "FP_adj_child"] = count[6]
            ti_p_df.loc[index, "FN_adj_child"] = count[7]

            # calculate delta rel. abundance - return [MSEsum, MSEsum_adj]
            MSEs = func.delta_relAbundance(sim_data=sim_ti_p, tool_data=row, children_data=json_data)

            ti_p_df.loc[index, "relAb_MSE"] = MSEs[0]
            ti_p_df.loc[index, "relAb_MSE_adj_child"] = MSEs[1]

        ti_p_df["timepoint"] = ti_p
        main_df = pd.concat([main_df, ti_p_df], ignore_index=True)

    meta_df = func.get_simulationStats(config, meta_df, exp_name, w_dir)
    main_df = main_df.merge(meta_df, on ="timepoint")
    main_df.reset_index(drop=True, inplace=True)

    print()
    print(main_df[['timepoint', 'tool_name', 
                    'TP', 
                    'TN', 
                    'FP', 
                    'FN', 
                    "relAb_MSE", 
                    "coverage_avg",
                    "coverage_sd",
                    "MAPQ_avg",
                    "uniformity_wg_per"
                ]])

    print()
    print(main_df[['timepoint', 'tool_name', 
                    "TP_adj_child",
                    "TN_adj_child", 
                    "FP_adj_child",
                    "FN_adj_child",
                    "relAb_MSE_adj_child",
                    "coverage_avg",
                    "coverage_sd",
                    "MAPQ_avg",
                    "uniformity_wg_per"
                ]])

    

    # export
    export_df = main_df[["timepoint", "tool_name", "sample",
                    "TP", "TN", "FP", "FN", 
                    "relAb_MSE", 
                    "TP_adj_child",
                    "TN_adj_child", 
                    "FP_adj_child",
                    "FN_adj_child",
                    "relAb_MSE_adj_child",
                    "coverage_avg",
                    "coverage_sd",
                    "MAPQ_avg",
                    "uniformity_wg_per"]].copy()
                    
    export_df.to_csv(w_dir / "experiments" / exp_name / "results" / "postPrediction" / "data.csv", index=False)


# execution starts here
parser = argparse.ArgumentParser()
parser.add_argument("-w", "--workding_directory", type=str, default='/home/asched85/mas_benchmark')
parser.add_argument("-e","--experiment", type=str, default="Ex01_WideQual Ex02_WideQual Ex03_WideQual Ex04_WideQual Ex05_WideQual Ex06_WideQual Ex07_WideQual Ex08_WideQual")
parser.add_argument("-p", "--parameter_file", type=str, default='workflow_config.yaml')
args=parser.parse_args()

w_dir = Path(args.workding_directory) 

with open(w_dir / args.parameter_file, 'r') as cfile:
    config = yaml.safe_load(cfile)

if len(args.experiment.split(" ")) == 1:
    exp = args.experiment.split(" ")[0].strip()
    print(f"Starting analysis for experiment {exp}")
    main(w_dir, config, exp_name=exp)
    
else:
    for exp in args.experiment.split(" "):
        print(f"Starting analysis for experiment {exp}")
        main(w_dir, config, exp_name=exp.strip())