from os.path import isdir
import re
import glob


def checkInput(config: dict) -> str:
    if config.get("input").lower().endswith("csv"):
        input_type = "trees"
    elif isdir(config.get("input")):
        folder = config.get("input")
        sheet = config.get("patientInfo", "NoSheetFound")
        files = glob.glob(f'{folder}/**/*', recursive=True)

        grn_idat_pattern = re.compile(r'Grn\.idat', re.IGNORECASE)
        red_idat_pattern = re.compile(r'Red\.idat', re.IGNORECASE)
        samplesheet_pattern = re.compile(rf'{sheet}', re.IGNORECASE)

        grn_idat_found = any(grn_idat_pattern.search(file) for file in files)
        red_idat_found = any(red_idat_pattern.search(file) for file in files)
        samplesheet_found = any(samplesheet_pattern.search(file) for file in files)

        if grn_idat_found and red_idat_found and samplesheet_found:
            input_type = "complete"
        elif grn_idat_found and red_idat_found and not samplesheet_found:
            raise RuntimeError("IDAT files detected, but the samplesheet couldn't be found in the input folder. Check that the file contains the 'samplesheet' keyword in its name and that it's a csv file.")
        else:
            raise RuntimeError(
                "Unexpected input folder. Check that the folder contains IDAT files for both red (Red.idat) and green (Grn.idat) and a samplesheet. Check that parameter --pattern matches your samplesheet. If the issue persists, contact the developers."
            )
    else:
        raise RuntimeError(
            """Couldn't determine the input type.
        - If running a single beta file, make sure the file is in csv format and ends with '.csv'
        - If running a folder of multiple csv files (multiple individuals), make sure csv files are in the folder 
        - If running the complete workflow, make sure a the samplesheet is in the folder (the file name should contain 'samplesheet') and that IDAT files are located recursively within the folder. 
        """
        )
    return input_type
