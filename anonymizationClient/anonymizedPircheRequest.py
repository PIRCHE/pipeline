import argparse
import requests
import csv
import json


def rest_request(url, username, password, api_key, api_req_payload):

    if not username or not password:
        payload = {'grant_type': 'authorization_code', 'client_id': 'api-client', 'code': api_key}
    else:
        payload = {'grant_type': 'password', 'client_id': 'portal-client', 'username': username, 'password': password}

    if verbose:
        print("auth_req_payload:", payload)

    r = requests.post(url + "/portal/rest/oauth/token", params=payload)

    if verbose:
        print('auth_url: ' + r.url)
        print('auth_req_status: ', r.status_code)
        print("auth_req_response: ", r.json())

    response = r.json()

    headers = {
        "Authorization": "Bearer " + response['access_token'],
        "Content-type": "application/json;charset=utf-8",
        "Accept": "application/json"
    }

    if verbose:
        print("api_req_payload:", json.dumps(api_req_payload))

    r = requests.post(url + "/pirche/rest/sot/api/match", headers=headers, json=api_req_payload)

    if verbose:
        print('api_url: ' + r.url)
        print('api_req_status: ', r.status_code)
        print("api_req_response: ", r.text)

    return r.text


def get_api_requests_payload(raw_input_data):

    donor_data = []
    api_requests = []

    for content in raw_input_data.values():
        tx_data = get_tx_data(content)
        patient_data = {'id': tx_data['patient']['id'], 'population': tx_data['patient']['population'], 'glString': tx_data['patient']['glString']}
        donor_data.append({'id': tx_data['donor']['id'], 'population': tx_data['donor']['population'], 'glString': tx_data['donor']['glString']})
        api_request = {'patient': patient_data, 'donors': donor_data}
        api_requests.append(api_request)
        donor_data = []

    return api_requests


def get_tx_data(raw_tx_data):

    p_id = raw_tx_data["pid"]
    d_id = raw_tx_data["did"]

    ploc_a = {'A1': raw_tx_data["pA1"], 'A2': raw_tx_data["pA2"]}
    ploc_b = {'B1': raw_tx_data["pB1"], 'B2': raw_tx_data["pB2"]}
    ploc_c = {'C1': raw_tx_data["pC1"], 'C2': raw_tx_data["pC2"]}
    ploc_drb = {'DRB11': raw_tx_data["pDRB11"], 'DRB12': raw_tx_data["pDRB12"]}
    ploc_dqb = {'DQB11': raw_tx_data["pDQB11"], 'DQB12': raw_tx_data["pDQB12"]}

    p_hla = {'A': ploc_a, 'B': ploc_b, 'C': ploc_c, 'DRB1': ploc_drb, 'DQB1': ploc_dqb}

    dloc_a = {'A1': raw_tx_data["dA1"], 'A2': raw_tx_data["dA2"]}
    dloc_b = {'B1': raw_tx_data["dB1"], 'B2': raw_tx_data["dB2"]}
    dloc_c = {'C1': raw_tx_data["dC1"], 'C2': raw_tx_data["dC2"]}
    dloc_drb = {'DRB11': raw_tx_data["dDRB11"], 'DRB12': raw_tx_data["dDRB12"]}
    dloc_dqb = {'DQB11': raw_tx_data["dDQB11"], 'DQB12': raw_tx_data["dDQB12"]}

    d_hla = {'A': dloc_a, 'B': dloc_b, 'C': dloc_c, 'DRB1': dloc_drb, 'DQB1': dloc_dqb}

    p_gl_string = build_gl_string(p_hla)
    d_gl_string = build_gl_string(d_hla)

    if not args.population:
        population = 'NMDP EUR haplotypes (2007)'
    else:
        population = args.population

    pat_data = {'id': p_id, 'population': population, 'hla': p_hla, 'glString': p_gl_string}
    don_data = {'id': d_id, 'population': population, 'hla': d_hla, 'glString': d_gl_string}
    tx_data = {'patient': pat_data, 'donor': don_data}

    return tx_data


def build_gl_string(hla):

    gl_string = ""
    for locus in hla.values():
        allele_count = 0
        for allele in locus.values():
            allele_count += 1
            if allele_count == 1:
                if allele:
                    gl_string += allele + "+"
            else:
                if allele:
                    gl_string += allele + "^"

    # remove trailing element
    gl_string = gl_string[:-1]

    return gl_string


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-url", "--url", help="target url", type=str, required=True)
    parser.add_argument("-u", "--user", help="username for request", type=str)
    parser.add_argument("-p", "--password", help="user password", type=str)
    parser.add_argument("-k", "--apikey", help="api key", type=str)
    parser.add_argument("-i", "--input", help="typing data input", type=str, required=True)
    parser.add_argument("-pp", "--population", help="population", type=str)

    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        print("verbose mode active")

    raw_csv_data = {}

    referenceHeader = 'pid,pA1,pA2,pB1,pB2,pC1,pC2,pDRB11,pDRB12,pDQB11,pDQB12,' \
                      'did,dA1,dA2,dB1,dB2,dC1,dC2,dDRB11,dDRB12,dDQB11,dDQB12'
    with open(args.input, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        rawHeader = csv_reader.fieldnames
        importHeader = res = ','.join([str(item) for item in rawHeader])
        if verbose:
            print(f'referenceHeader: {referenceHeader}')
            print(f'importHeader: {importHeader}')
        if referenceHeader == importHeader:
            line_count = 0
            for row in csv_reader:
                if verbose:
                    print(f'raw_csv_row: {row["pid"]} {row["pA1"]} {row["pA2"]} {row["pB1"]} {row["pB2"]} {row["pC1"]}'
                          f' {row["pC2"]} {row["pDRB11"]} {row["pDRB12"]} {row["pDQB11"]} {row["pDQB12"]}'
                          f' {row["did"]} {row["dA1"]} {row["dA2"]} {row["dB1"]} {row["dB2"]} {row["dC1"]}'
                          f' {row["dC2"]} {row["dDRB11"]} {row["dDRB12"]} {row["dDQB11"]} {row["dDQB12"]}')
                raw_csv_data[line_count] = row
                line_count += 1
            print(f'Processed {line_count} lines.')
        else:
            print("CSV file provided does not match expected format.")

    api_requests_payload = get_api_requests_payload(raw_csv_data)

    for api_request_payload in api_requests_payload:
        rest_request(args.url, args.user, args.password, args.apikey, api_request_payload)
