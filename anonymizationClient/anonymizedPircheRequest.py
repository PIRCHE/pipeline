import argparse
import requests
import csv
import json
import xlrd
import random
import operator


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

    return r


def get_api_requests_payload(raw_input_data, genotype_data):

    api_donor_data = []
    api_requests = []

    for content in raw_input_data.values():
        tx_data = get_tx_data(content)
        patient_data = {'id': tx_data['patient']['id'], 'population': tx_data['patient']['population'], 'glString': tx_data['patient']['glString']}
        api_donor_data.append(({'id': tx_data['donor']['id'], 'population': tx_data['donor']['population'], 'glString': tx_data['donor']['glString']}))
        if args.anonymization:
            fk_data = get_fk_data(tx_data, genotype_data)
            api_requests.extend(fk_data)
        api_request = {'patient': patient_data, 'donors': api_donor_data}
        api_requests.append(api_request)
        api_donor_data = []

    return api_requests


def get_tx_data(raw_tx_data):

    check_unsupported_loci(raw_tx_data)

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


def get_fk_data(tx_data, genotype_data):

    fk_data = []
    num_fk_genotypes = 5

    p_hla_tx = tx_data['patient']['hla']
    d_hla_tx = tx_data['donor']['hla']

    genotypes = genotype_data['genotypes']

    p_genotypes = get_genotype_subset(p_hla_tx, genotypes, num_fk_genotypes)
    d_genotypes = get_genotype_subset(d_hla_tx, genotypes, num_fk_genotypes)

    for idx in range(num_fk_genotypes):
        api_don_data = []
        pat_data = genotype_to_person_data(p_genotypes[idx])
        api_don_data.append(genotype_to_person_data(d_genotypes[idx]))
        fk_data.append({'patient': pat_data, 'donors': api_don_data})

    return fk_data


def get_genotype_subset(hla_tx, genotypes, k_genotypes):

    hlas = []

    for loci in hla_tx.values():
        # TODO handle/ignore empty loci/allele value
        hlas.append(loci)

    hlas_sel = random.sample(hlas, 2)

    for idx, item in enumerate(hlas_sel[0].items()):
        if idx == 0:
            ref_loc1_k1 = item[0]
            ref_loc1_v1 = item[1]
        elif idx == 1:
            ref_loc1_k2 = item[0]
            ref_loc1_v2 = item[1]

    for idx, item in enumerate(hlas_sel[1].items()):
        if idx == 0:
            ref_loc2_k1 = item[0]
            ref_loc2_v1 = item[1]
        elif idx == 1:
            ref_loc2_k2 = item[0]
            ref_loc2_v2 = item[1]

    if verbose:
        print('refloc1_key1: ' + ref_loc1_k1 + ' -- refloc1_value1: ' + ref_loc1_v1)
        #print('refloc1_key2: ' + ref_loc1_k2 + ' -- refloc1_value2: ' + ref_loc1_v2)
        #print('refloc2_key1: ' + ref_loc2_k1 + ' -- refloc2_value1: ' + ref_loc2_v1)
        #print('refloc2_key2: ' + ref_loc2_k2 + ' -- refloc2_value2: ' + ref_loc2_v2)

        print('genotype_len: ', len(genotypes))

    genotype_subset = [loci for loci in genotypes if loci[ref_loc1_k1] == ref_loc1_v1]
                       #and loci[ref_loc1_k2] == ref_loc1_v2]
                       #and loci[ref_loc2_k1] == ref_loc2_v1]
                       #and loci[ref_loc2_k2] == ref_loc2_v2]

    if verbose:
        print('genotype_subset_len: ', len(genotype_subset))
        for idx, genotype in enumerate(genotype_subset):
            if idx == 15:
                break
            print('subset genotypes: ', genotype)

    #genotype_subset.sort(key=operator.itemgetter('freq'), reverse=True)

    # print('after sort: ', genotype_subset)
    # if verbose:
    #     for idx, genotype in enumerate(genotype_subset):
    #         if idx == 15:
    #             break
    #         print('after sort: ', genotype)

    frequencies = [genotype['freq'] for genotype in genotype_subset]
    if verbose:
        for idx, frequency in enumerate(frequencies):
            if idx == 15:
                break
            print('subset frequencies: ', frequency)

    genotypes_random = random.choices(genotype_subset, frequencies, k=k_genotypes)

    if verbose:
        print(genotypes_random)

    return genotypes_random


def genotype_to_person_data(genotype):

    id = genotype["id"]

    loc_a = {'A1': genotype["A1"], 'A2': genotype["A2"]}
    loc_b = {'B1': genotype["B1"], 'B2': genotype["B2"]}
    loc_c = {'C1': genotype["C1"], 'C2': genotype["C2"]}
    loc_drb = {'DRB11': genotype["DRB11"], 'DRB12': genotype["DRB12"]}
    loc_dqb = {'DQB11': genotype["DQB11"], 'DQB12': genotype["DQB12"]}

    hla = {'A': loc_a, 'B': loc_b, 'C': loc_c, 'DRB1': loc_drb, 'DQB1': loc_dqb}

    print('fake_gl_string_input:', hla)

    gl_string = build_gl_string(hla)

    if not args.population:
        population = 'NMDP EUR haplotypes (2007)'
    else:
        population = args.population

    person_data = {'id': id, 'population': population, 'hla': hla, 'glString': gl_string}

    return person_data

def check_unsupported_loci(raw_tx_data):
    if raw_tx_data["pDRB31"] or raw_tx_data["pDRB32"] \
            or raw_tx_data["pDRB41"] or raw_tx_data["pDRB42"] \
            or raw_tx_data["pDRB51"] or raw_tx_data["pDRB52"] \
            or raw_tx_data["pDQA11"] or raw_tx_data["pDQA12"] \
            or raw_tx_data["pDPA11"] or raw_tx_data["pDPA12"] \
            or raw_tx_data["pDPB11"] or raw_tx_data["pDPB12"] \
            or raw_tx_data["dDRB31"] or raw_tx_data["dDRB32"] \
            or raw_tx_data["dDRB41"] or raw_tx_data["dDRB42"] \
            or raw_tx_data["dDRB51"] or raw_tx_data["dDRB52"] \
            or raw_tx_data["dDQA11"] or raw_tx_data["dDQA12"] \
            or raw_tx_data["dDPA11"] or raw_tx_data["dDPA12"] \
            or raw_tx_data["dDPB11"] or raw_tx_data["dDPB12"]:
        print('INFO: Only the locus A, B, C, DRB1 and DQB1 are supported '
              'the hla data of other locus will be ignored.')


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


def cleanHLA(locus, allele):
    allele = allele.replace("g", "")
    if "*" in allele:
        allele = allele[allele.index("*") + 1:]
    if ":" not in allele:
        if len(allele) > 3:
            allele = allele[0:2] + ":" + allele[2:4]
        else:
            print("ERROR: cannot clean " + locus + " " + allele)
    return locus + "*" + allele.replace("g", "")


def build_genotypes():

    # TODO handle g groups
    # read files
    spacer = '-'*40
    if verbose:
        print("read haplotype file")
        print(spacer)

    workbook = xlrd.open_workbook(args.haplotypes)
    worksheet = workbook.sheet_by_index(0)

    header = worksheet.row(0)
    freq_suffix = "_freq"
    rank_suffix = "_rank"

    col_idx = {}
    for idx, cell_obj in enumerate(header):
        if verbose:
            print(str(idx) + " -> " + cell_obj.value)
        col_idx[cell_obj.value] = idx
    if verbose:
        print(spacer)

    haplotypes = []

    for row_idx in range(1, worksheet.nrows):
        if worksheet.cell(row_idx, col_idx[args.population_short + rank_suffix]).value != "NA":
            haplotype = {
                'id': row_idx,
                'A':  cleanHLA("A", worksheet.cell(row_idx, col_idx["A"]).value),
                'B':  cleanHLA("B", worksheet.cell(row_idx, col_idx["B"]).value),
                'C':  cleanHLA("C", worksheet.cell(row_idx, col_idx["C"]).value),
                'DRB1':  cleanHLA("DRB1", worksheet.cell(row_idx, col_idx["DRB1"]).value),
                'DQB1':  cleanHLA("DQB1", worksheet.cell(row_idx, col_idx["DQB1"]).value),
                'freq':  float(worksheet.cell(row_idx, col_idx[args.population_short + freq_suffix]).value),
            }
            haplotypes.append(haplotype)

    if verbose:
        print("read " + str(len(haplotypes)) + " haplotypes")
        print(spacer)

    haplotypes.sort(key=lambda x: x['freq'], reverse=True)

    genotypes = []
    frequencies = []
    genotypes_count = 0
    threshold = float(args.threshold)
    cumulated_left = 0.
    #with open(args.output, "w") as output:
        #writer = csv.writer(output, delimiter=',')
        #writer.writerow(["genotype", "A1", "A2", "B1", "B2", "C1", "C2", "DRB11", "DRB12", "DQB11", "DQB12", "freq"])
    for left in haplotypes:
        cumulated_left += left['freq']
        #if cumulated_left < threshold and verbose:
            #print("Remaining frequency " + "{:.5f}".format(threshold - cumulated_left))
        cumulated_right = 0.
        for right in haplotypes:
            cumulated_right += right['freq']
            if cumulated_left < threshold and cumulated_right < threshold and left['id'] != right['id']:
                #writer.writerow([str(left['id']) + "-" + str(right['id']), left['A'], right['A'], left['B'], right['B'], left['C'], right['C'], left['DRB1'], right['DRB1'], left['DQB1'], right['DQB1'], left['freq'] * right['freq']])
                genotype = {'id': str(left['id']) + "-" + str(right['id']), 'A1': left['A'], 'A2': right['A'], 'B1': left['B'], 'B2': right['B'], 'C1': left['C'], 'C2': right['C'], 'DRB11': left['DRB1'], 'DRB12': right['DRB1'], 'DQB11': left['DQB1'], 'DQB12': right['DQB1'], 'freq': left['freq'] * right['freq']}
                genotypes.append(genotype)
                genotypes_count += 1
                frequencies.append(left['freq'] * right['freq'])

    if verbose:
        print(spacer)
        print("created " + str(genotypes_count) + " genotypes based on threshold of " + str(threshold))
        print(spacer)
        print("successfully stored genotypes in file")
        print(spacer)

    #frequencies.sort(key=operator.itemgetter('freq'), reverse=True)
    genotypes.sort(key=operator.itemgetter('freq'), reverse=True)

    # gt_count = 0
    # for gt in genotypes:
    #     print(gt)
    #     gt_count += 1
    #     if gt_count == 150:
    #         break

    genotype_data = {'genotypes': genotypes, 'frequencies': frequencies}

    return genotype_data


def write_results(match_results):
    with open(args.output, "w") as output:
        writer = csv.writer(output, delimiter=',')
        writer.writerow(["id", "p1_A", "p1_B", "p1_C", "p1_DRB1", "p1_DRB3", "p1_DRB4", "p1_DRB5", "p1_DPA1", "p1_DPB1", "p1_DQA1", "p1_DQB1", "p1_SUM"
                            , "p2_A", "p2_B", "p2_C", "p2_DRB1", "p2_DRB3", "p2_DRB4", "p2_DRB5", "p2_DPA1", "p2_DPB1", "p2_DQA1", "p2_DQB1", "p2_SUM"])

        for match_result in match_results:
            p1_scores = match_result['pircheI_scores']
            p2_scores = match_result['pircheII_scores']
            writer.writerow([str(match_result['id']), p1_scores['A'], p1_scores['B'], p1_scores['C'], p1_scores['DRB1'], p1_scores['DRB3'], p1_scores['DRB4'], p1_scores['DRB5'], p1_scores['DPA1'], p1_scores['DPB1'], p1_scores['DQA1'], p1_scores['DQB1'], p1_scores['sum']
                            , p2_scores['A'], p2_scores['B'], p2_scores['C'], p2_scores['DRB1'], p2_scores['DRB3'], p2_scores['DRB4'], p2_scores['DRB5'], p2_scores['DPA1'], p2_scores['DPB1'], p2_scores['DQA1'], p2_scores['DQB1'], p2_scores['sum']])


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-url", "--url", help="target url", type=str, required=True)
    parser.add_argument("-u", "--user", help="username for request", type=str)
    parser.add_argument("-p", "--password", help="user password", type=str)
    parser.add_argument("-k", "--apikey", help="api key", type=str)
    parser.add_argument("-i", "--input", help="typing data input", type=str, required=True)
    parser.add_argument("-pp", "--population", help="population", type=str)
    parser.add_argument("-hp", "--haplotypes", help="NMDP haplotype table (either 2007 or 2011 or equally formatted)", required=True)
    parser.add_argument("-t", "--threshold", help="frequency threshold for haplotypes 0.0 to 1.0", required=True)
    parser.add_argument("-ps", "--population_short", help="population short code as used in the header row", required=True)
    parser.add_argument("-o", "--output", help="output file name", required=True)
    parser.add_argument("-a", "--anonymization", help="Enable anonymization. Default - False - no anonymization", default=False)

    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        print("verbose mode active")

    if args.anonymization:
        genotype_data = build_genotypes()
    else:
        genotype_data = []

    raw_csv_data = {}

    referenceHeader = 'pid,pA1,pA2,pB1,pB2,pC1,pC2,pDRB11,pDRB12,pDRB31,pDRB32,pDRB41,pDRB42,pDRB51,pDRB52,' \
                      'pDQA11,pDQA12,pDQB11,pDQB12,pDPA11,pDPA12,pDPB11,pDPB12,' \
                      'did,dA1,dA2,dB1,dB2,dC1,dC2,dDRB11,dDRB12,dDRB31,dDRB32,dDRB41,dDRB42,dDRB51,dDRB52,' \
                      'dDQA11,dDQA12,dDQB11,dDQB12,dDPA11,dDPA12,dDPB11,dDPB12'

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

    api_requests_payload = get_api_requests_payload(raw_csv_data, genotype_data)

    results = []
    for api_request_payload in api_requests_payload:
        response = rest_request(args.url, args.user, args.password, args.apikey, api_request_payload)
        response_raw = response.json()
        response_raw_p1 = response_raw["pircheI"]
        response_raw_p2 = response_raw["pircheII"]
        result_data = {'id': list(response_raw_p1.keys())[0], 'pircheI_scores': list(response_raw_p1.values())[0], 'pircheII_scores': list(response_raw_p2.values())[0]}
        results.append(result_data)

    write_results(results)
