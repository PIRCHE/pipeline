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

    input_count = 0
    donor_data = []
    api_requests = []
    for content in raw_input_data.values():
        if (input_count % 2) == 0:
            patient_data = get_person_data(content)
        else:
            donor_data.append(get_person_data(content))
            api_request = {'patient': patient_data, 'donors': donor_data}
            api_requests.append(api_request)
            donor_data = []
        input_count += 1

    return api_requests


def get_person_data(raw_person_data):

    row_count = 0
    tuple_count = 1
    gl_string = ""
    for row_element in raw_person_data.values():
        if row_count == 0:
            person_id = row_element
        else:
            if tuple_count == 1:
                if row_element:
                    gl_string += row_element + "+"
                tuple_count += 1
            else:
                if row_element:
                    gl_string += row_element + "^"
                tuple_count = 1
        row_count += 1

    # remove trailing element
    gl_string = gl_string[:-1]

    if not args.population:
        population = 'NMDP EUR haplotypes (2007)'
    else:
        population = args.population

    person_data = {'id': person_id, 'population': population, 'glString': gl_string}

    return person_data


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

    referenceHeader = 'id,A1,A2,B1,B2,C1,C2,DRB11,DRB12,DQB11,DQB12'
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
                    print(f'raw_csv_row: {row["id"]} {row["A1"]} {row["A2"]} {row["B1"]} {row["B2"]} {row["C1"]}'
                          f' {row["C2"]} {row["DRB11"]} {row["DRB12"]} {row["DQB11"]} {row["DQB12"]}')
                raw_csv_data[line_count] = row
                line_count += 1
            print(f'Processed {line_count} lines.')
        else:
            print("CSV file provided does not match expected format.")

    api_requests_payload = get_api_requests_payload(raw_csv_data)

    for api_request_payload in api_requests_payload:
        rest_request(args.url, args.user, args.password, args.apikey, api_request_payload)
