import argparse
import requests
import csv
import json
import xlrd
import random
import operator
import hashlib
import uuid
import itertools

# TODO
# low resolution support


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


def get_api_requests_data(raw_input_data, genotype_data):

    api_donor_data = []
    api_requests_data = []

    for content in raw_input_data.values():
        tx_data = get_tx_data(content)
        id_map = {'p_id': (tx_data['patient']['id'], str(uuid.uuid4())), 'd_id': (tx_data['donor']['id'], str(uuid.uuid4()))}
        patient_data = {'id': id_map['p_id'][1], 'population': tx_data['patient']['population'], 'glString': tx_data['patient']['glString']}
        api_donor_data.append(({'id': id_map['d_id'][1], 'population': tx_data['donor']['population'], 'glString': tx_data['donor']['glString']}))
        if args.anonymization:
            fk_data = get_fk_data(tx_data, genotype_data)
            for fk_data_entry in fk_data:
                api_request_data = {'api_payload': fk_data_entry}
                api_requests_data.append(api_request_data)
        api_request_payload = {'patient': patient_data, 'donors': api_donor_data}
        api_request_data = {'api_payload': api_request_payload, 'id_map': id_map}
        api_requests_data.append(api_request_data)
        random.shuffle(api_requests_data)
        api_donor_data = []

    return api_requests_data


def get_tx_data(raw_tx_data):

    check_unsupported_loci(raw_tx_data)

    p_id = raw_tx_data["pid"]
    d_id = raw_tx_data["did"]

    ploc_a = {'A1': next(iter(clean_hla("A", raw_tx_data["pA1"]))), 'A2': next(iter(clean_hla("A", raw_tx_data["pA2"])))}
    ploc_b = {'B1': next(iter(clean_hla("B", raw_tx_data["pB1"]))), 'B2': next(iter(clean_hla("B", raw_tx_data["pB2"])))}
    ploc_c = {'C1': next(iter(clean_hla("C", raw_tx_data["pC1"]))), 'C2': next(iter(clean_hla("C", raw_tx_data["pC2"])))}
    ploc_drb = {'DRB11': next(iter(clean_hla("DRB1", raw_tx_data["pDRB11"]))), 'DRB12': next(iter(clean_hla("DRB1", raw_tx_data["pDRB12"])))}
    ploc_dqb = {'DQB11': next(iter(clean_hla("DQB1", raw_tx_data["pDQB11"]))), 'DQB12': next(iter(clean_hla("DQB1", raw_tx_data["pDQB12"])))}

    p_hla = {'A': ploc_a, 'B': ploc_b, 'C': ploc_c, 'DRB1': ploc_drb, 'DQB1': ploc_dqb}

    dloc_a = {'A1': next(iter(clean_hla("A", raw_tx_data["dA1"]))), 'A2': next(iter(clean_hla("A", raw_tx_data["dA2"])))}
    dloc_b = {'B1': next(iter(clean_hla("B", raw_tx_data["dB1"]))), 'B2': next(iter(clean_hla("B", raw_tx_data["dB2"])))}
    dloc_c = {'C1': next(iter(clean_hla("C", raw_tx_data["dC1"]))), 'C2': next(iter(clean_hla("C", raw_tx_data["dC2"])))}
    dloc_drb = {'DRB11': next(iter(clean_hla("DRB1", raw_tx_data["dDRB11"]))), 'DRB12': next(iter(clean_hla("DRB1", raw_tx_data["dDRB12"])))}
    dloc_dqb = {'DQB11': next(iter(clean_hla("DQB1", raw_tx_data["dDQB11"]))), 'DQB12': next(iter(clean_hla("DQB1", raw_tx_data["dDQB12"])))}

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
    num_fk_genotypes = int(args.kanonymization)

    p_glString_tx = tx_data['patient']['glString']
    d_glString_tx = tx_data['donor']['glString']

    genotypes = genotype_data['genotypes']

    p_genotypes = convert_genotypes(get_fk_genotypes(p_glString_tx, genotypes, num_fk_genotypes, "pat"))
    d_genotypes = convert_genotypes(get_fk_genotypes(d_glString_tx, genotypes, num_fk_genotypes, "don"))

    p_genotypes.append(tx_data['patient']['hla'])
    d_genotypes.append(tx_data['donor']['hla'])
    gt_matrix = itertools.product(p_genotypes, d_genotypes)

    for gt in gt_matrix:
        if gt[0] != tx_data['patient']['hla'] or gt[1] != tx_data['donor']['hla']:  # exclude real patient donor pair here
            api_don_data = []
            pat_data = genotype_to_person_data(gt[0])
            api_don_data.append(genotype_to_person_data(gt[1]))
            fk_data.append({'patient': pat_data, 'donors': api_don_data})

    return fk_data


def get_fk_genotypes(hla_tx, genotypes, k_genotypes, person_type):

    genotypes_random = []

    k_gen_real = random.randrange(1, k_genotypes-1)
    k_gen_fk = k_genotypes - k_gen_real

    salt = args.salt
    salted_input = hla_tx + salt + person_type

    md5_hex = hashlib.md5(salted_input.encode()).hexdigest()
    md5_str = str(md5_hex)

    frequencies = [genotype['freq'] for genotype in genotypes]

    # if verbose:
    #     for idx, frequency in enumerate(frequencies):
    #         if idx == 15:
    #             break
    #         print('frequencies: ', frequency)

    random.seed(md5_str)
    genotypes_real = random.choices(genotypes, frequencies, k=k_gen_real)
    genotypes_random.extend(genotypes_real)

    genotypes_fk = build_fake_genotypes(genotypes, frequencies, k_gen_fk)
    genotypes_random.extend(genotypes_fk)

    if verbose:
        print(genotypes_random)

    return genotypes_random


def build_fake_genotypes(genotypes, frequencies, k_gen_fk):

    fk_genotypes = []
    genotypes_fk_rnd = random.choices(genotypes, frequencies, k=k_gen_fk*5)

    for i in range(k_gen_fk):
        gen_fk = {'id': '',
                  'A1': genotypes_fk_rnd[i]["A1"], 'A2': genotypes_fk_rnd[i]["A2"],
                  'B1': genotypes_fk_rnd[i*2]["B1"], 'B2': genotypes_fk_rnd[i*2]["B2"],
                  'C1': genotypes_fk_rnd[i*3]["C1"], 'C2': genotypes_fk_rnd[i*3]["C2"],
                  'DRB11': genotypes_fk_rnd[i*4]["DRB11"], 'DRB12': genotypes_fk_rnd[i*4]["DRB12"],
                  'DQB11': genotypes_fk_rnd[i*5]["DQB11"], 'DQB12': genotypes_fk_rnd[i*5]["DQB12"],
                  'freq': genotypes_fk_rnd[i]["freq"]}
        fk_genotypes.append(gen_fk)

    return fk_genotypes


def convert_genotypes(genotypes):

    genotypes_converted = []
    for genotype in genotypes:
        loc_a = {'A1': next(iter(genotype["A1"])), 'A2': next(iter(genotype["A2"]))}
        loc_b = {'B1': next(iter(genotype["B1"])), 'B2': next(iter(genotype["B2"]))}
        loc_c = {'C1': next(iter(genotype["C1"])), 'C2': next(iter(genotype["C2"]))}
        loc_drb = {'DRB11': next(iter(genotype["DRB11"])), 'DRB12': next(iter(genotype["DRB12"]))}
        loc_dqb = {'DQB11': next(iter(genotype["DQB11"])), 'DQB12': next(iter(genotype["DQB12"]))}
        genotypes_converted.append({'A': loc_a, 'B': loc_b, 'C': loc_c, 'DRB1': loc_drb, 'DQB1': loc_dqb})

    return genotypes_converted


def genotype_to_person_data(genotype):

    p_id = str(uuid.uuid4())

    if verbose:
        print('fake_gl_string_input:', genotype)

    gl_string = build_gl_string(genotype)

    if not args.population:
        population = 'NMDP EUR haplotypes (2007)'
    else:
        population = args.population

    person_data = {'id': p_id, 'population': population, 'glString': gl_string}

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
    for locus, locus_full in hla.items():
        allele_count = 0
        # print('locus_full:', locus_full)
        for allele in locus_full.values():
            allele_count += 1
            if allele_count == 1:
                if allele:
                    if "*" in allele:
                        gl_string += allele + "+"
                    else:
                        gl_string += locus + "*" + allele + "+"
                    # print('gl_string_build:', gl_string)
            else:
                if allele:
                    if "*" in allele:
                        gl_string += allele + "^"
                    else:
                        gl_string += locus + "*" + allele + "^"
                    # print('gl_string_build:', gl_string)

    # remove trailing element
    gl_string = gl_string[:-1]

    return gl_string


def clean_hla(locus, allele):
    alleles = []
    if "g" in allele:
        # alleles = g_groups.get(locus+allele)  # 2005 file format
        alleles = g_groups.get(allele)  # 2011 file format
    else:
        if "*" in allele:
            allele = allele[allele.index("*") + 1:]
        if ":" not in allele:
            if len(allele) > 3:
                allele = allele[0:2] + ":" + allele[2:4]
            else:
                print("ERROR: cannot clean " + locus + " " + allele)
        alleles.append(allele)
    return alleles


def build_genotypes():

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
                'A':  clean_hla("A", worksheet.cell(row_idx, col_idx["A"]).value),
                'B':  clean_hla("B", worksheet.cell(row_idx, col_idx["B"]).value),
                'C':  clean_hla("C", worksheet.cell(row_idx, col_idx["C"]).value),
                'DRB1':  clean_hla("DRB1", worksheet.cell(row_idx, col_idx["DRB1"]).value),
                'DQB1':  clean_hla("DQB1", worksheet.cell(row_idx, col_idx["DQB1"]).value),
                'freq':  float(worksheet.cell(row_idx, col_idx[args.population_short + freq_suffix]).value),
            }
            haplotypes.append(haplotype)

    if verbose:
        print("read " + str(len(haplotypes)) + " haplotypes")
        print(spacer)

    haplotypes.sort(key=lambda x: x['freq'], reverse=True)

    genotypes = []
    genotypes_count = 0
    threshold = float(args.threshold)
    cumulated_left = 0.
    for left in haplotypes:
        cumulated_left += left['freq']
        cumulated_right = 0.
        for right in haplotypes:
            cumulated_right += right['freq']
            if cumulated_left < threshold and cumulated_right < threshold and left['id'] != right['id']:
                genotype = {'id': str(left['id']) + "-" + str(right['id']), 'A1': left['A'], 'A2': right['A'], 'B1': left['B'], 'B2': right['B'], 'C1': left['C'], 'C2': right['C'], 'DRB11': left['DRB1'], 'DRB12': right['DRB1'], 'DQB11': left['DQB1'], 'DQB12': right['DQB1'], 'freq': left['freq'] * right['freq']}
                genotypes.append(genotype)
                genotypes_count += 1

    if verbose:
        print(spacer)
        print("created " + str(genotypes_count) + " genotypes based on threshold of " + str(threshold))
        print(spacer)
        print("successfully stored genotypes in file")
        print(spacer)

    genotypes.sort(key=operator.itemgetter('freq'), reverse=True)

    # if verbose:
    #     gt_count = 0
    #     for gt in genotypes:
    #         print(gt)
    #         gt_count += 1
    #         if gt_count == 150:
    #             break

    genotype_data = {'genotypes': genotypes}

    return genotype_data


def build_g_groups():
    g_groups = {}
    with open(args.ggroups) as f:
        for line in f:
            line_strip = line.rstrip()
            line_values = line_strip.split(';')
            if len(line_values) > 2 and line_values[2] != '':
                group_alleles_unsort = clean_g_group_alleles(line_values[1])
                group_alleles_sort = list(group_alleles_unsort)
                group_alleles_sort.sort()
                # group_name = "".join(line_values[2].split(":", 2)[:2]) + 'g'  # 2005 haplo file format
                # locus = line_values[0].replace('*', '')  # 2005 haplo file format
                group_name = ":".join(line_values[2].split(":", 2)[:2]) + 'g'  # 2011 haplo file format
                locus = line_values[0]  # 2011 haplo file format
                g_groups[locus + group_name] = group_alleles_sort
                print('locus:', locus, ' -- group_alleles:', group_alleles_sort, ' -- group_name:', group_name)
    return g_groups


def clean_g_group_alleles(g_group_alleles):
    group_alleles = set()
    alleles = g_group_alleles.split('/')
    for allele in alleles:
        if 'N' not in allele and 'Q' not in allele:
            clean_allele = ":".join(allele.split(":", 2)[:2])
            if clean_allele not in group_alleles:
                group_alleles.add(clean_allele)
    return group_alleles


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
    parser.add_argument("-gg", "--ggroups", help="HLA g-groups reference table", required=True)
    parser.add_argument("-hp", "--haplotypes", help="NMDP haplotype table (either 2007 or 2011 or equally formatted)", required=True)
    parser.add_argument("-t", "--threshold", help="frequency threshold for haplotypes 0.0 to 1.0", required=True)
    parser.add_argument("-ps", "--population_short", help="population short code as used in the header row", required=True)
    parser.add_argument("-o", "--output", help="output file name", required=True)
    parser.add_argument("-a", "--anonymization", help="Enable anonymization. Default - False - no anonymization", default=False)
    parser.add_argument("-s", "--salt", help="Salt (password) used to anonymize input data. Use identical password when submitting same HLA data set multiple times.")
    parser.add_argument("-ka", "--kanonymization", help="Number of smoke hla data sets (genotypes) generated per patient and per donor. Default - 5.", default=5)

    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        print("verbose mode active")

    g_groups = build_g_groups()

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

    api_requests_data = get_api_requests_data(raw_csv_data, genotype_data)

    results = []
    for api_request_data in api_requests_data:
        response = rest_request(args.url, args.user, args.password, args.apikey, api_request_data["api_payload"])
        response_raw = response.json()
        response_raw_p1 = response_raw["pircheI"]
        response_raw_p2 = response_raw["pircheII"]
        if "id_map" in api_request_data and len(response_raw_p1.keys()) > 0 and len(response_raw_p1.values()) and len(response_raw_p2.values()) > 0:
            if list(response_raw_p1.keys())[0] == api_request_data["id_map"]["d_id"][1]:
                result_data = {'id': api_request_data["id_map"]["d_id"][0], 'pircheI_scores': list(response_raw_p1.values())[0], 'pircheII_scores': list(response_raw_p2.values())[0]}
                results.append(result_data)
            else:
                print('ERROR: request(' + api_request_data["id_map"]["d_id"][1] + ') and response (' + list(response_raw_p1.keys())[0] + ') ids do not match. Result for donor_ID (' + api_request_data["id_map"]["d_id"][0] + ') skipped.')

    write_results(results)
