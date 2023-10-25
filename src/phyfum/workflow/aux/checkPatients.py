import pandas as pd


def checkPatientsBetas(patients: list[str], dir: str) -> list[str]:
    validPatients = []
    for patient in patients:
        csv_file = f"{dir}/{patient}.csv"

        df = pd.read_csv(csv_file)
        if len(df.columns) >= 3:
            validPatients.append(patient)
    if len(validPatients) == 0:
        raise ValueError("No patients with at least 3 samples were found! Please check your input data.")

    return validPatients


def checkPatientInfo(patientInfo, patient_col):
    patientInfo = patientInfo[patientInfo.groupby(patient_col)[patient_col].transform('count') >= 3]
    if len(patientInfo) == 0:
        raise ValueError("No patients with at least 3 samples were found! Please check your input data.")
