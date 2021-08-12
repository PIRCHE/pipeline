import argparse
import requests
import json

def rest_request(url, username):
    r =requests.get(url + "/blabla")
    response = json.loads(r.text)

    headers = {
        "Authorization": "Bearer " + response.access_token
    }
    r =requests.get(url + "/authorized_request", headers = headers)

    return r.text


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-i", "--input", help="typing data input", type=str, required=True)
    parser.add_argument("-u", "--user", help="username for request", type=str, required=True)
    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        print("verbose mode active")

    rest_request("https://www.pirche.com", args.user)


    input_data = {}

    with open(args.input, 'r') as input_file:
        for row in input_file:
            cols = row.split(",")
            print(cols)
            input_data[cols[0]] = cols