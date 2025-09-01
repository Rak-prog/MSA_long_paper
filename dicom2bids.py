import pydicom
import nibabel as nb
import argparse 
import numpy as np
from loguru import logger
import pydicom.tag
from utils.utils import load_img, save_img, parse_yaml, str_to_bool, LoguruStreamHandler # this is the package with the functions shared across project which you can find and edit here: /home/riccardo/codici_progetti/shared_utils
import sys, os, yaml, re
from collections import Counter, defaultdict
from pathlib import Path
import csv
import pandas as pd
import subprocess 

def calc_time_between_timepoints(df):
    rows = []

    # Group by subject
    for subj_id, group in df.groupby("Subj_ID"):
        group = group.sort_values("Session").reset_index(drop=True)

        # Extract session -> date mapping
        session_dates = dict(zip(group["Session"], group["ExamDate"]))
        
        # Compute time differences
        t0 = session_dates.get("ses-T0")
        t1 = session_dates.get("ses-T1")
        t2 = session_dates.get("ses-T2")

        for i, row in group.iterrows():
            diff_t1_t0 = (t1 - t0).days if t1 and t0 else None
            diff_t2_t1 = (t2 - t1).days if t2 and t1 else None
            diff_t2_t0 = (t2 - t0).days if t2 and t0 else None

            rows.append({
                "Subj_ID": row["Subj_ID"],
                "Session": row["Session"],
                "ExamDate": row["ExamDate"],
                "T1_T0_Diff_Days": diff_t1_t0,
                "T2_T1_Diff_Days": diff_t2_t1,
                "T2_T0_Diff_Days": diff_t2_t0
            })

    # Build new dataframe
    print(rows)
    df_diff = pd.DataFrame(rows)
    return df_diff

def is_ordered(sessions):
    expected_order = ['ses-T0', 'ses-T1', 'ses-T2','ses-T3']
    return list(sessions) == sorted(sessions, key=lambda x: expected_order.index(x))

# Group by Subj_ID and Session and check if dates are strictly increasing
def check_monotonic(group):
    return group["ExamDate"].is_monotonic_increasing

def parse_args(args):
    """parse argument"""
    parser = argparse.ArgumentParser(epilog = "Sript to convert the dicom data to BIDS, running BIDS validation and if dataset is valid runs also optionally mriqc.\
                                      It prepares also an excel with all the dicom tags information. and it checks the date in case of longitudinal studies")
    parser.add_argument("-i", "--input", required=True, type = str, help = 'The input must be a file with the following columns T1_path, Subj ID, Session')
    parser.add_argument("-i_clinical","--input_clinical", required = False, type = str, help = 'This is the excel file wit the clinical information, such as disease duration and clinical scales')
    parser.add_argument("-o", "--output", required=True, type = str, help = "OutputFolder where to create the BIDS root")
    parser.add_argument("-of", "--output_file", required=True, type = str, help = "Name of output file in .xlsx format")
    parser.add_argument("-rc", "--runconversion", required=True, type = str_to_bool, default = True, help = "True if you want to run the dicom to bids conversion. If True also output an excel file with all the dicom tag extracted (slower). \
                                                                                                    So think carefully iy you want to deactivate this. I usually set this to False when I want only to run BIDS validator")
    parser.add_argument("-dp", "--doplot", required=True, type = str_to_bool, default = False, help = "True/False, if True visualize the just converted nifti")
    parser.add_argument("-dv", "--dovalidator", required=True, type = str_to_bool, default = False, help = "True/False, if True runs teh docker for validating if the new dataset is BIDS compliant")
    parser.add_argument("-dqc", "--doqualitycheck", required=True, type = str_to_bool, default = False, help = "True/False, if True runs the docker for MRIQC. This can be runned only if validotr retunrs a green flag")
    parser.add_argument("-list", "--subjectlist", required=True, nargs = "+", help = "IF quality check is True you must pass here a list of subjects onto which perform the quality check, you can avoid passing the prefix sub-")
    parser.add_argument("-group", "--groupmriqc", required=True, type = str_to_bool, default = False, help = "True/False, if True runs the docker for MRIQC FOR GROUP LEVEL ANALYSIS AND REPORTS. This can be runned only if validotr retunrs a green flag and the sinle level has been peformed")

    return parser.parse_args()

def main(argv):
    try:
        args = parse_args(argv)
        #config = parse_yaml(args)
        if "parameters_check_SP_LAST_CORRETTA_TUTTO" in args.input:
            data = pd.read_excel(args.input,sheet_name="parameters_check",index_col=None, engine='openpyxl')
        elif "controlli" in args.input: # if you have the files to convert the controls
            data = pd.read_excel(args.input,index_col=None, engine='openpyxl')

        if args.input_clinical and "parameters_check_SP_LAST_CORRETTA_TUTTO" in args.input:  # do emrge with clinical data only if the argument is passed in the command line and for the group of patients
            print(args.input_clinical)
            data_clinical = pd.read_excel(args.input_clinical, sheet_name = 'pazienti_MSA_all', index_col = False, engine = 'openpyxl')
            #data.to_csv('prova2.csv', index = False)
            # Merge only ses-T0 rows
            t0_df = data[data['Session'] == 'ses-T0'].merge(data_clinical, on='Subj_ID', how='left')
            non_t0_df = data[data['Session'] != 'ses-T0']
            data = pd.concat([t0_df, non_t0_df], ignore_index=True)
            data = data.sort_values(['Subj_ID', 'Session']).reset_index(drop=True)
            #data.to_csv('prova.csv', index = False)
            

        result = data.groupby('Subj_ID')['Session'].apply(is_ordered)
        print(result)
        quality_check_output=f"{args.output}/derivatives/mriqc"
        quality_check_workdir=f"{args.output}/derivatives/mriqc_workdir"
        print(data.head(3))

        if args.runconversion == True:
            # new columns 
            data["ExamDate"] = None
            data["FieldStrength"] = None
            data["RepetitionTime"] = None
            data["EchoTime"] = None
            data["FlipAngle"] = None
            data["MRAcquisitionType"] = None
            data["Series"] = None
            data["Model"] = None 
            data["Manufacturer"] = None 
            data["SeriesDescription"] = None
            data["Slices"] = None
            data["Thickness"] = None
            data["SpacingX"] = None
            data["SpacingY"] = None
            data["Age"] = None
            data["Sex"] = None
            data["Height"] = None
            data["Weight"] = None
            data["Orientation"] = None  # after conversion get neurological or radiological
            data["InstitutionName"] = None
            data["InstitutionAddress"] = None
            data["StationName"] = None
            data["DeviceSerialNumber"] = None
            data["SoftwareVersions"] = None
            
            for index, row in data.iterrows():
                t1_path = row["T1_path"]
                #subj_name = row["Subj_Name"]
                subj_id = row["Subj_ID"]
                session = row["Session"]

                try:
                    dicom_files = [f for f in os.listdir(t1_path) if not f.startswith(".")]

                    # anonimizza qua tutti i dicom presenti nella cartella 
                    
                    dcm = pydicom.dcmread(os.path.join(t1_path, dicom_files[0]), stop_before_pixels=True)
                    #print(dcm)
                    if not dicom_files:
                        print(f"No files in {t1_path}")
                        continue
                    data.at[index, "ExamDate"] = getattr(dcm, "StudyDate", None)
                    data.at[index, "FieldStrength"] = getattr(dcm, "MagneticFieldStrength", None)
                    data.at[index, "RepetitionTime"] = getattr(dcm, "RepetitionTime", None)
                    data.at[index, "EchoTime"] = getattr(dcm, "EchoTime", None)
                    data.at[index, "FlipAngle"] = getattr(dcm, "FlipAngle", None)
                    data.at[index, "MRAcquisitionType"] = getattr(dcm, "MRAcquisitionType", None)
                    data.at[index, "Series"] = getattr(dcm, "SeriesDescription", None)
                    data.at[index, "Model"] = getattr(dcm, "ManufacturerModelName", None) 
                    data.at[index, "Manufacturer"] = getattr(dcm, "Manufacturer", None) 
                    data.at[index, "SeriesDescription"] = getattr(dcm, "SeriesDescription", None) 
                    data.at[index, "Slices"] = int(len(dicom_files)) 
                    data.at[index, "Thickness"] = getattr(dcm, "SliceThickness", None) 
                    data.at[index, "SpacingX"] = getattr(dcm, "PixelSpacing", None)[0]
                    data.at[index, "SpacingY"] = getattr(dcm, "PixelSpacing", None)[1]
                    data.at[index,"Age"] = getattr(dcm, "PatientAge", None) 
                    data.at[index,"Sex"] = getattr(dcm, "PatientSex", None) 
                    data.at[index,"Height"] = getattr(dcm, "PatientSize", None) 
                    data.at[index,"Weight"] = getattr(dcm, "PatientWeight", None) 
                    data.at[index,"InstitutionName"] = getattr(dcm, "InstitutionName", None) 
                    data.at[index,"InstitutionAddress"] = getattr(dcm, "InstitutionAddress", None) 
                    data.at[index,"StationName"] = getattr(dcm,"StationName",None)
                    data.at[index,"DeviceSerialNumber"] = getattr(dcm, "DeviceSerialNumber", None) 
                    data.at[index,"SoftwareVersions"] = getattr(dcm,"SoftwareVersions",None)
                    data.at[index,"ParallelReductionFactorInPlane"] = getattr(dcm,"ParallelReductionFactorInPlane",None)
                    data.at[index,"MatrixCoilMode"] = getattr(dcm,"MatrixCoilMode",None)
                    data.at[index,"ReceiveCoilName"] = getattr(dcm,"ReceiveCoilName",None)

                except Exception as e:
                    print(f"Error reading DICOM in {t1_path}: {e}")
                
                # create the root tree            
                subject_dir=Path(f"{args.output}/{subj_id}/{session}/anat")
                if not os.path.exists(subject_dir):
                    os.makedirs(subject_dir)
    
                #print(f"{args.output}/{subj_id}/{session}")
                #print(t1_path)

                # to avoid re-running it given that it goes throug a loop
                if os.path.isfile(f"{subject_dir}/{subj_id}_{session}_T1w.nii.gz"):
                    print(f"subject {subj_id} and session {session} already converted")
                else:
                    subprocess.run(["dcm2niix", "-m", "y", "-z", "y", "-ba", "y", "-o", str(subject_dir), "-f", f"{subj_id}_{session}_T1w" , t1_path])       # subprocess.run

                if subj_id == "sub-18" and session == "ses-T0": # this is the only subjects acquired with a Pihilips 1.5 T it is also in PSR orientation not RAS 
                    print("doing sub-18, ses-T0")
                    if os.path.isfile(f"{subject_dir}/{subj_id}_{session}_T1w_real.nii.gz"):
                        subprocess.run(["fslreorient2std", f"{subject_dir}/{subj_id}_{session}_T1w_real.nii.gz", f"{subject_dir}/{subj_id}_{session}_T1w.nii.gz"])
                        os.remove(f"{subject_dir}/{subj_id}_{session}_T1w_real.nii.gz")
                        os.rename(f"{subject_dir}/{subj_id}_{session}_T1w_real.json",f"{subject_dir}/{subj_id}_{session}_T1w.json")
                        #subprocess.run(["fslswapdim", f"{str(subject_dir)}/{subj_id}_{session}_T1w.nii.gz","LR", "AP", "IS", f"{str(subject_dir)}/{subj_id}_{session}_T1w.nii.gz"])
                        subprocess.run(["fsleyes", f"{str(subject_dir)}/{subj_id}_{session}_T1w.nii.gz"])

                hdr , img = load_img(f"{subject_dir}/{subj_id}_{session}_T1w.nii.gz")    
                data.at[index,"Orientation"] = nb.aff2axcodes(hdr.affine)

                if args.doplot == True:
                    subprocess.run(["fsleyes", f"{str(subject_dir)}/{subj_id}_{session}_T1w.nii.gz"])

                print("#########")

            #data.to_excel("MSA_new_Dataset.xlsx", index=False)

            # check date for longitudinal
            data["ExamDate"] = pd.to_datetime(data["ExamDate"].astype(str), format="%Y%m%d")
            df_sorted = data.sort_values(by=["Subj_ID", "Session", "ExamDate"])
            check =( df_sorted.groupby(["Subj_ID", "Session"])["ExamDate"].apply(lambda x: x.is_monotonic_increasing))
            #print(check)
            non_monotonic = check[check == False]

            if non_monotonic.empty:
                print("All dates are in ascending order per subject and session.")
            else:
                print("Some subject/session groups have non-ascending dates:")
                #print(non_monotonic)
            
            # calc_time_between_timepoints
            diff = calc_time_between_timepoints(data.sort_values(by=["Subj_ID", "Session"]))
            #print(diff)

            print(diff.shape)
            print(data.shape)

            merged = pd.merge(data, diff, on=["Subj_ID","Session","ExamDate"], how="left") # suffixes=("", "_dup"))
            
            merged['years'] = merged.apply(lambda row: 0 if row['Session'] == 'ses-T0' else row['T1_T0_Diff_Days'], axis=1)
            merged.to_excel(args.output_file, index=False)

            #printing some interesting statistics to the terminal 
            for statistics in ['Session','RepetitionTime','EchoTime','Series','Slices']:
                print(f"Grouping by {statistics} ################################")
                session_counts = merged.groupby('Subj_ID')[statistics].nunique()
                count_distribution = session_counts.value_counts().sort_index()
                print(count_distribution)
                print("################################")
            
            # check if ses-T0 and ses-T1 are in the correct order
            result = merged.groupby('Subj_ID')['Session'].apply(is_ordered)
            print(result)

        if args.dovalidator == True:
            subprocess.run(["sudo", "docker", "run", "-ti", "--rm", "-v", f"{args.output}:/data:ro", "bids/validator"  ,"/data" ,"--ignoreWarnings", "--ignoreSymlinks" ,"--ignoreNiftiHeaders", "--ignoreSubjectConsistency", "--verbose"]) #--verbose
        if args.doqualitycheck == True:
            list_pass = " ".join(args.subjectlist)
            if args.groupmriqc == False:
                subprocess.run(["sudo", "docker", "run", "-ti", "--rm", "-v", f"{args.output}:/data:ro", "-w", f"{quality_check_workdir}:/work:ro", "-v", f"{quality_check_output}:/out",  \
                                "nipreps/mriqc:latest", "/data", "/out", "participant",\
                                "--n_cpus","16","--omp-nthreads", "8","--mem_gb","32","--participant-label",f"{list_pass}","--write-graph","--verbose-reports","--no-sub"]) #--verbose
            elif args.groupmriqc == True:
                subprocess.run(["sudo", "docker", "run", "-ti", "--rm", "-v", f"{args.output}:/data:ro", "-w", f"{quality_check_workdir}:/work:ro", "-v", f"{quality_check_output}:/out",  \
                                "nipreps/mriqc:latest", "/data", "/out", "group","--participant-label",f"{list_pass}"])

    except Exception as e:
        logger.error(f"An error occurred during the script: {e}")
    finally:
        logger.info("Finished processing data.")

if __name__ == "__main__":
    sys.exit(main(sys.argv))

