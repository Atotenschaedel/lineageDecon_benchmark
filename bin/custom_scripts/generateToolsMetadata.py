# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 22:30:22 2023
This script generates meta data file for vaquero tool based on parameter from config file
@author: aschedl
"""
# import libraries
import os, yaml
import pandas as pd
import argparse

# # add flag argument
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str)
parser.add_argument("-t","--tool_name", type=str)
parser.add_argument("-o","--output", type=str)
args=parser.parse_args()


def vaquero(df, config: dict):
    """
    Function to generate metadata file for vaquero input

    Parameters
    ----------
    df   :  pandas dataframe containing "sample" and "sample_date" column
    config: dictonary containing additional metadata 

    Returns
    -------
    df: pandas dataframe
    """
    df["BSF_run"] = config["VAQUERO"]["BSF_RUN"]
    df["BSF_start_date"] = config["VAQUERO"]["BSF_START_DATE"]
    df["LocationID"] = config["VAQUERO"]["LOCATION_ID"]
    df["LocationName"] = config["VAQUERO"]["LOCATION_NAME"]
    df["N_in_Consensus"] = config["VAQUERO"]["N_IN_CONSENSUS"]
    df["adress_town"] = config["VAQUERO"]["ADDRESS_TOWN"]
    df["connected_people"] = config["VAQUERO"]["CONNECTED_PEOPLE"]
    df["dcpLatitude"] = config["VAQUERO"]["DCP_LATITUDE"]
    df["dcpLongitude"] = config["VAQUERO"]["DCP_LONGTITUTDE"]
    df["include_in_report"] = config["VAQUERO"]["INCLUDE_IN_REPORT"]
    df["report_category"] = config["VAQUERO"]["REPORT_CATEGORY"]
    df["status"] = config["VAQUERO"]["STATUS"]
    df["additional_information"] = ""
    df["BSF_sample_name"] = df["sample"] + config["VAQUERO"]["SAMPLE_SUFFIX"]

    df.rename(columns={"sample":"RNA_ID_int"}, inplace=True)
    return df

# get parameter file
with open(os.getcwd()+"/workflow_config.yaml", "r") as cfile:
    config = yaml.safe_load(cfile)

df = pd.read_csv(args.input, sep="\t")

# choose tool to generate data
if args.tool_name == "vaquero":
    df_export=vaquero(df, config)

df_export.to_csv(args.output, index=False, sep="\t")