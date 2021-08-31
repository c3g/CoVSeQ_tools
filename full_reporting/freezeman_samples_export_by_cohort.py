import argparse
import csv
import subprocess
import os
import requests

# This script is a rudimentary example to fetch data from the Freezeman API.
# This extract a limited amount of fields from the API for a single cohort, but more are available.
# (check on the Freezeman GUI API page to browse the endpoints)
# Launch this script from the command line :
# $ python samples_export_by_cohort.py -url http://biobank.genome.mcgill.ca/api/ -user FMS_USER -password FMS_PASSWORD -csv OUTPUT_FILE -cohort FILTERING_PARAM

#FMS_API_BASE_URL = "http://biobank.genome.mcgill.ca/api/"
# cohort HostSeq = BQC
# cohort VirusSeq = INSPQ_COVID

AUTH_TOKEN_ENDPOINT = "token/"
SAMPLES_ENDPOINT = "samples/list_export/"

def execute(fms_base_url, fms_user, fms_password, output_file, fms_cert, cohort):
    # Set up proxies to access internet : uncomment the 2 following lines to run the script from inside the center.
    #os.environ["http_proxy"] = "http://192.168.32.1:3128"
    #os.environ["https_proxy"] = "http://192.168.32.1:3128"
    
    # Setup Certificate CA_BUNDLE location
    PATH_TO_CERT = fms_cert

    print("Starting extraction of sample data from freezeman")

    # Get jwt token
    print("Requesting authorization...")
    auth = requests.post(fms_base_url + AUTH_TOKEN_ENDPOINT, data={"username": fms_user, "password": fms_password}, verify=PATH_TO_CERT)
    if auth.status_code == 200:
        access_token = auth.json()["access"]
        headers = {"Authorization": "Bearer " + access_token}

        # Extract data
        print("Requesting data...")
        data = []

        params = {"individual__cohort": cohort, "format": "csv"}
        response = requests.get(args.url + SAMPLES_ENDPOINT, params=params, headers=headers, verify=PATH_TO_CERT)

        if response.status_code == 200:
            print("Receiving data...")
            try:
                with open(output_file, mode="w") as f:
                    f.write(response.content.decode('utf-8'))
            except Exception as e:
                print("Failed writing result file : " + str(e.message) + " (" + e.args + ")")

            print("Operation complete!")
        else:
            print("Failed call to Freezeman API : " + str(response.status_code) + " (" + response.text + ")")
    else:
        print("Failed to authenticate...")
        print(auth.text)

if __name__ == '__main__':
    # Get parameters from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-csv", help="Path of the output csv file")
    parser.add_argument("-url", help="Freezeman API base url")
    parser.add_argument("-user", help="Freezeman User")
    parser.add_argument("-password", help="Freezeman Password")
    parser.add_argument("-cert", help="Location of CA Bundle certification")
    parser.add_argument("-cohort", help="Cohort of the individual from which the sample has been extracted")
    args = parser.parse_args()

    execute(args.url, args.user, args.password, args.csv, args.cert, args.cohort)
