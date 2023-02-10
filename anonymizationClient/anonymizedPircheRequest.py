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


def rest_request(url, username, password, api_key, api_req_payload, proxies):

    headers = auth_request(url, username, password, api_key, proxies)

    if verbose:
        print("api_req_payload:", json.dumps(api_req_payload))

    if args.requestmode == "multi":
        r = requests.post(url + "/pirche/rest/sot/multipatientmatch", headers=headers, proxies=proxies, data=api_req_payload)
    elif args.requestmode == "single":
        r = requests.post(url + "/pirche/rest/sot/api/match", headers=headers, proxies=proxies, json=api_req_payload)

    if verbose:
        print('api_url: ' + r.url)
        print('api_req_status: ', r.status_code)
        print("api_req_response: ", r.text)

    return r


def request_population_id(url, username, password, api_key, proxies, p_name):

    headers = auth_request(url, username, password, api_key, proxies)

    pop_r = requests.get(url + "/pirche/rest/populations", headers=headers, proxies=proxies)

    for population in pop_r.json():
        if population["name"] == p_name:
            p_id = population["id"]

    if verbose:
        print('api_url: ' + pop_r.url)
        print('api_req_status: ', pop_r.status_code)
        print("api_req_response: ", pop_r.text)

    return p_id


def auth_request(url, username, password, api_key, proxies):

    if dev:
        print('proxies: ', proxies)

    if not username or not password:
        payload = {'grant_type': 'authorization_code', 'client_id': 'api-client', 'code': api_key}
    else:
        payload = {'grant_type': 'password', 'client_id': 'portal-client', 'username': username, 'password': password}

    if verbose:
        print("auth_req_payload:", payload)

    auth_r = requests.post(url + "/portal/rest/oauth/token", proxies=proxies, params=payload)

    if verbose:
        print('auth_url: ' + auth_r.url)
        print('auth_req_status: ', auth_r.status_code)
        print("auth_req_response: ", auth_r.json())

    auth_response = auth_r.json()

    headers = {
        "Authorization": "Bearer " + auth_response['access_token'],
        "Content-type": "application/json;charset=utf-8",
        "Accept": "application/json"
    }

    return headers

def get_api_request_data_multi(raw_input_data, genotype_data, proxies):

    api_request_data = {}
    api_payloads = []
    id_list = {}
    full_data = []

    count_datasets = 0

    pop_id = request_population_id(args.url, args.user, args.password, args.apikey, proxies, args.population)

    header = "0.5,500,1000,0.01,0," + pop_id + ",false,false,false,false,false,false,2,2,2,2,2,150,true,true,false,0,0,false|"
    body = ""

    for content in raw_input_data.values():

        tx_data = get_tx_data(content)
        id_map = {'p_id': (tx_data['patient']['id'], str(uuid.uuid4())), 'd_id': (tx_data['donor']['id'], str(uuid.uuid4()))}
        id_list[id_map['d_id'][1]] = id_map
        patient_data_tx = {'id': id_map['p_id'][1], 'population': tx_data['patient']['population'], 'hlaString': tx_data['patient']['hlaString']}
        donor_data_tx = {'id': id_map['d_id'][1], 'population': tx_data['donor']['population'], 'hlaString': tx_data['donor']['hlaString']}
        tx_pair = {"patient": patient_data_tx, "donor": donor_data_tx}
        full_data.append(tx_pair)

        if args.anonymization:
            fk_data = get_fk_data(tx_data, genotype_data)
            full_data.extend(fk_data)

    random.shuffle(full_data)

    for data_entry in full_data:
        body += data_entry["patient"]["id"] + "," + data_entry["patient"]["hlaString"] + "\n"
        body += data_entry["donor"]["id"] + "," + data_entry["donor"]["hlaString"] + "\n"
        body += "," + "\n"

        count_datasets += 1

        # if dev:
        #     print('Count:', count_datasets)

        if count_datasets == args.requestsize:
            api_payload = header + body.rsplit("\n", 2)[0]
            api_payloads.append(api_payload)
            body = ""
            count_datasets = 0

    if count_datasets > 0:
        api_payload = header + body.rsplit("\n", 2)[0]
        api_payloads.append(api_payload)

    api_request_data = {'api_payloads': api_payloads, 'id_list': id_list}

    return api_request_data

def get_api_requests_data_single(raw_input_data, genotype_data):
    api_donor_data = []
    api_requests_data = []

    for content in raw_input_data.values():
        tx_data = get_tx_data(content)
        id_map = {'p_id': (tx_data['patient']['id'], str(uuid.uuid4())), 'd_id': (tx_data['donor']['id'], str(uuid.uuid4()))}
        patient_data = {'id': id_map['p_id'][1], 'population': tx_data['patient']['population'], 'glString': tx_data['patient']['hlaString']}
        api_donor_data.append(({'id': id_map['d_id'][1], 'population': tx_data['donor']['population'], 'glString': tx_data['donor']['hlaString']}))
        if args.anonymization:
            fk_data = get_fk_data(tx_data, genotype_data)
            for fk_data_entry in fk_data:
                api_request_data = {'api_payload': fk_data_entry}  #create on multi csv request content object
                api_requests_data.append(api_request_data)
        api_request_payload = {'patient': patient_data, 'donors': api_donor_data}
        api_request_data = {'api_payload': api_request_payload, 'id_map': id_map}
        api_requests_data.append(api_request_data)
        api_donor_data = []

    random.shuffle(api_requests_data)

    if dev:
        print('number datasets:', len(api_requests_data))

    return api_requests_data

def get_tx_data(raw_tx_data):

    check_unsupported_loci(raw_tx_data)
    check_nmdp_codes(raw_tx_data)

    p_id = raw_tx_data["pid"]
    d_id = raw_tx_data["did"]

    ploc_a = {'alleles': {'A1': next(iter(clean_hla("A", raw_tx_data["pA1"]))), 'A2': next(iter(clean_hla("A", raw_tx_data["pA2"])))}, 'res': check_loc_res([raw_tx_data["pA1"], raw_tx_data["pA2"]])}
    ploc_b = {'alleles': {'B1': next(iter(clean_hla("B", raw_tx_data["pB1"]))), 'B2': next(iter(clean_hla("B", raw_tx_data["pB2"])))}, 'res': check_loc_res([raw_tx_data["pB1"], raw_tx_data["pB2"]])}
    ploc_c = {'alleles': {'C1': next(iter(clean_hla("C", raw_tx_data["pC1"]))), 'C2': next(iter(clean_hla("C", raw_tx_data["pC2"])))}, 'res': check_loc_res([raw_tx_data["pC1"], raw_tx_data["pC2"]])}
    ploc_drb = {'alleles': {'DRB11': next(iter(clean_hla("DRB1", raw_tx_data["pDRB11"]))), 'DRB12': next(iter(clean_hla("DRB1", raw_tx_data["pDRB12"])))}, 'res': check_loc_res([raw_tx_data["pDRB11"], raw_tx_data["pDRB12"]])}
    ploc_dqb = {'alleles': {'DQB11': next(iter(clean_hla("DQB1", raw_tx_data["pDQB11"]))), 'DQB12': next(iter(clean_hla("DQB1", raw_tx_data["pDQB12"])))}, 'res': check_loc_res([raw_tx_data["pDQB11"], raw_tx_data["pDQB12"]])}

    p_hla = {'A': ploc_a, 'B': ploc_b, 'C': ploc_c, 'DRB1': ploc_drb, 'DQB1': ploc_dqb}

    dloc_a = {'alleles': {'A1': next(iter(clean_hla("A", raw_tx_data["dA1"]))), 'A2': next(iter(clean_hla("A", raw_tx_data["dA2"])))}, 'res': check_loc_res([raw_tx_data["dA1"], raw_tx_data["dA2"]])}
    dloc_b = {'alleles': {'B1': next(iter(clean_hla("B", raw_tx_data["dB1"]))), 'B2': next(iter(clean_hla("B", raw_tx_data["dB2"])))}, 'res': check_loc_res([raw_tx_data["dB1"], raw_tx_data["dB2"]])}
    dloc_c = {'alleles': {'C1': next(iter(clean_hla("C", raw_tx_data["dC1"]))), 'C2': next(iter(clean_hla("C", raw_tx_data["dC2"])))}, 'res': check_loc_res([raw_tx_data["dC1"], raw_tx_data["dC2"]])}
    dloc_drb = {'alleles': {'DRB11': next(iter(clean_hla("DRB1", raw_tx_data["dDRB11"]))), 'DRB12': next(iter(clean_hla("DRB1", raw_tx_data["dDRB12"])))}, 'res': check_loc_res([raw_tx_data["dDRB11"], raw_tx_data["dDRB12"]])}
    dloc_dqb = {'alleles': {'DQB11': next(iter(clean_hla("DQB1", raw_tx_data["dDQB11"]))), 'DQB12': next(iter(clean_hla("DQB1", raw_tx_data["dDQB12"])))}, 'res': check_loc_res([raw_tx_data["dDQB11"], raw_tx_data["dDQB12"]])}

    d_hla = {'A': dloc_a, 'B': dloc_b, 'C': dloc_c, 'DRB1': dloc_drb, 'DQB1': dloc_dqb}

    if args.requestmode == "multi":
        p_hla_string = build_csv_data(p_hla)
        d_hla_string = build_csv_data(d_hla)
    else:
        p_hla_string = build_gl_string(p_hla)
        d_hla_string = build_gl_string(d_hla)

    if not args.population:
        population = 'NMDP EUR haplotypes (2007)'
    else:
        population = args.population

    pat_data = {'id': p_id, 'population': population, 'hla': p_hla, 'hlaString': p_hla_string}
    don_data = {'id': d_id, 'population': population, 'hla': d_hla, 'hlaString': d_hla_string}
    tx_data = {'patient': pat_data, 'donor': don_data}

    return tx_data


def get_fk_data(tx_data, genotype_data):

    fk_data = []
    num_fk_genotypes = int(args.kanonymization)

    p_hlaString_tx = tx_data['patient']['hlaString']
    d_hlaString_tx = tx_data['donor']['hlaString']

    genotypes = genotype_data['genotypes']

    p_genotypes = convert_genotypes(get_fk_genotypes(p_hlaString_tx, genotypes, num_fk_genotypes, "pat"))
    d_genotypes = convert_genotypes(get_fk_genotypes(d_hlaString_tx, genotypes, num_fk_genotypes, "don"))

    handle_res(tx_data['patient']['hla'], p_genotypes)
    handle_res(tx_data['donor']['hla'], d_genotypes)

    p_genotypes.append(tx_data['patient']['hla'])
    d_genotypes.append(tx_data['donor']['hla'])
    gt_matrix = itertools.product(p_genotypes, d_genotypes)

    for gt in gt_matrix:
        if gt[0] != tx_data['patient']['hla'] or gt[1] != tx_data['donor']['hla']:  # exclude real patient donor pair here
            pat_data = genotype_to_person_data(gt[0])
            if args.requestmode == "multi":
                don_data = genotype_to_person_data(gt[1])
                fk_data.append({'patient': pat_data, 'donor': don_data})
            elif args.requestmode == "single":
                don_data = [genotype_to_person_data(gt[1])]
                fk_data.append({'patient': pat_data, 'donors': don_data})

        else:
            if dev:
                print('real tx pair excluded pat: ', gt[0])
                print('real tx pair excluded donor: ', gt[1])

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

        if random.random() < 0.01:
            locus = random.choice(list(rare_alleles.keys()))
            allele = random.choice(list(rare_alleles[locus]))
            key = locus + random.choice(["1", "2"])
            gen_fk[key] = allele
        fk_genotypes.append(gen_fk)
    return fk_genotypes


def convert_genotypes(genotypes):

    genotypes_converted = []
    for genotype in genotypes:
        loc_a = {'alleles': {'A1': next(iter(genotype["A1"])), 'A2': next(iter(genotype["A2"]))}}
        loc_b = {'alleles': {'B1': next(iter(genotype["B1"])), 'B2': next(iter(genotype["B2"]))}}
        loc_c = {'alleles': {'C1': next(iter(genotype["C1"])), 'C2': next(iter(genotype["C2"]))}}
        loc_drb = {'alleles': {'DRB11': next(iter(genotype["DRB11"])), 'DRB12': next(iter(genotype["DRB12"]))}}
        loc_dqb = {'alleles': {'DQB11': next(iter(genotype["DQB11"])), 'DQB12': next(iter(genotype["DQB12"]))}}
        genotypes_converted.append({'A': loc_a, 'B': loc_b, 'C': loc_c, 'DRB1': loc_drb, 'DQB1': loc_dqb})

    return genotypes_converted


def genotype_to_person_data(genotype):

    p_id = str(uuid.uuid4())

    if dev:
        print('fake genotypes build data_string input:', genotype)

    hlaStringIdentifier = "hlaString"
    if args.requestmode == "multi":
        hla_string = build_csv_data(genotype)
    elif args.requestmode == "single":
        hla_string = build_gl_string(genotype)
        hlaStringIdentifier = "glString"

    if dev:
        print('fake genotypes build hlaString result:', hla_string)

    if not args.population:
        population = 'NMDP EUR haplotypes (2007)'
    else:
        population = args.population

    person_data = {'id': p_id, 'population': population, hlaStringIdentifier: hla_string}

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
        print('INFO: Only the locus A, B, C, DRB1 and DQB1 are supported. '
              'The hla data of other locus will be ignored.')

def check_nmdp_codes(raw_tx_data):
    if is_nmdp_code(raw_tx_data["pA1"]) or is_nmdp_code(raw_tx_data["pA2"]) \
            or is_nmdp_code(raw_tx_data["pB1"]) or is_nmdp_code(raw_tx_data["pB2"]) \
            or is_nmdp_code(raw_tx_data["pC1"]) or is_nmdp_code(raw_tx_data["pC2"]) \
            or is_nmdp_code(raw_tx_data["pDRB11"]) or is_nmdp_code(raw_tx_data["pDRB12"]) \
            or is_nmdp_code(raw_tx_data["pDQB11"]) or is_nmdp_code(raw_tx_data["pDQB12"]):
        print('ERROR: HLA typings containing NMDP codes are not supported.')
        exit()

def is_nmdp_code(raw_tx_data_value):
    is_nmdp = False
    if ":" in raw_tx_data_value:
        check_split = raw_tx_data_value.split(":", 1)
        if check_split[1].isalpha():
            is_nmdp = True
    return is_nmdp

def build_gl_string(hla):

    gl_string = ""
    for locus, locus_full in hla.items():
        allele_count = 0
        for allele in locus_full["alleles"].values():
            allele_count += 1
            if allele_count == 1:
                if allele:
                    if "*" in allele:
                        gl_string += allele + "+"
                    else:
                        gl_string += locus + "*" + allele + "+"
            else:
                if allele:
                    if "*" in allele:
                        gl_string += allele + "^"
                    else:
                        gl_string += locus + "*" + allele + "^"

    # remove trailing element
    gl_string = gl_string[:-1]

    return gl_string

def build_csv_data(hla):
    csv_data = ""

    for locus, locus_full in hla.items():
        for allele in locus_full["alleles"].values():
            if allele:
                csv_data += locus + "*" + allele + ","

    # remove trailing element
    csv_data = csv_data[:-1]

    return csv_data


def clean_hla(locus, allele):
    alleles = []
    if "g" in allele:
        if "*" in allele:  # e.g. g-group values format in 2011 haplotype file
            alleles = g_groups_global.get(allele)
        if "*" not in allele:  # e.g. g-group values format in 2007 haplotype file
            alleles = g_groups_global.get(locus+allele)
        if dev:
            print('g-group alleles', alleles)
    else:
        if "*" in allele:
            allele = allele[allele.index("*") + 1:]

        loc_prefixes = ["A", "B", "C", "DRB", "DQ", "DP"]
        if "*" not in allele:
            # order in list is important otherwise B will be removed from DRB before causing DR which is not removed afterwards
            pref_to_remove = {'DRB1': '', 'DRB3': '', 'DRB4': '', 'DRB5': '', 'DQA1': '', 'DQB1': '', 'DPA1': '', 'DPB1': '', 'A': '', 'B': '', 'C': ''}
            if any(loc_prefix in allele for loc_prefix in loc_prefixes):
                for key, value in pref_to_remove.items():
                    allele = allele.replace(key, value)
            else:
                if ":" not in allele:  # no '*' and no ':' e.g. allele format in 2007 haplotype file
                    if len(allele) > 3:
                        allele = allele[0:2] + ":" + allele[2:4]

        alleles.append(allele)
    return alleles


def check_loc_res(allele_values):
    loc_res = "mol_high"
    alleles_not_empty = [allele for allele in allele_values if allele.strip()]
    if any(":" not in allele for allele in alleles_not_empty):
        if any("*" not in allele for allele in alleles_not_empty):
            loc_res = 'ser'
        else:
            loc_res = 'mol_low'
    return loc_res


def handle_res(tx_hla_data, fk_gts):

    locs_res = []

    for loc, locus_data in tx_hla_data.items():
        locs_res.append({"loc": loc, "res": locus_data["res"]})
        if dev:
            print('handleRes - loc: ', loc)
            print('handleRes - locus_data: ', locus_data)

    for fk_gt in fk_gts:
        for loc_res in locs_res:
            fk_gt_loc = fk_gt[loc_res["loc"]]
            fk_gt_loc["res"] = loc_res["res"];
            if dev:
                print('fk_gt_loc', fk_gt_loc)
                print('loc_res', loc_res)
            if "ser" in loc_res["res"]:
                alleles_ser = {}
                for allele_loc, allele_val in fk_gt_loc["alleles"].items():
                    if dev:
                        print('allele to map to ser', allele_val)
                    if str(loc_res["loc"]) + allele_val in dna_ser_table_global.keys():
                        alleles_ser[allele_loc] = dna_ser_table_global[str(loc_res["loc"]) + allele_val]
                    else: # fallback to molecular low if no serological equivalent could be found due to mapping backwards from high to low
                        if dev:
                            print('fallback used to map to ser no serological equivalent found for', str(loc_res["loc"]) + allele_val)
                        alleles_ser[allele_loc] = (str(allele_val).split(":", 1)[0])
                fk_gt_loc["alleles"] = alleles_ser
                if dev:
                    print('ser mapped alleles', fk_gt_loc["alleles"])
            elif "mol_low" in loc_res["res"]:
                alleles_mol_low = {}
                for allele_loc, allele_val in fk_gt_loc["alleles"].items():
                    if dev:
                        print('allele to map to low', allele_val)
                    alleles_mol_low[allele_loc] = (str(allele_val).split(":", 1)[0])
                fk_gt_loc["alleles"] = alleles_mol_low
                if dev:
                    print('low mapped alleles', fk_gt_loc["alleles"])

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
        if dev:
            print("haplo file header: " + str(idx) + " -> " + cell_obj.value)
        col_idx[cell_obj.value] = idx
    if verbose:
        print(spacer)

    haplotypes = []

    for row_idx in range(1, worksheet.nrows):
        if worksheet.cell(row_idx, col_idx[args.haplofilepop + rank_suffix]).value != "NA":
            haplotype = {
                'id': row_idx,
                'A':  clean_hla("A", worksheet.cell(row_idx, col_idx["A"]).value),
                'B':  clean_hla("B", worksheet.cell(row_idx, col_idx["B"]).value),
                'C':  clean_hla("C", worksheet.cell(row_idx, col_idx["C"]).value),
                'DRB1':  clean_hla("DRB1", worksheet.cell(row_idx, col_idx["DRB1"]).value),
                'DQB1':  clean_hla("DQB1", worksheet.cell(row_idx, col_idx["DQB1"]).value),
                'freq':  float(worksheet.cell(row_idx, col_idx[args.haplofilepop + freq_suffix]).value),
            }
            haplotypes.append(haplotype)

    if verbose:
        print("read " + str(len(haplotypes)) + " haplotypes")
        print(spacer)

    haplotypes.sort(key=lambda x: x['freq'], reverse=True)

    genotypes = []
    genotypes_count = 0
    threshold = float(args.haplothreshold)
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

    # if dev:
    #     gt_count = 0
    #     for gt in genotypes:
    #         print(gt)
    #         gt_count += 1
    #         if gt_count == 150:
    #             break

    genotype_data = {'genotypes': genotypes}

    return genotype_data


def build_g_groups():
    g_grps = {}
    with open(args.ggroups) as f:
        for line in f:
            if not line.startswith("#"):
                line_strip = line.rstrip()
                line_values = line_strip.split(';')
                if len(line_values) > 2 and line_values[2] != '':
                    group_alleles_unsort = clean_g_group_alleles(line_values[1])
                    group_alleles_sort = list(group_alleles_unsort)
                    group_alleles_sort.sort()
                    if args.haplofileformat == "XXXX":  # 2007 haplo file format
                        group_name = "".join(line_values[2].split(":", 2)[:2]) + 'g'
                        locus = line_values[0].replace('*', '')
                    elif args.haplofileformat == "L*XX:XX":  # 2011 haplo file format
                        group_name = ":".join(line_values[2].split(":", 2)[:2]) + 'g'
                        locus = line_values[0]
                    g_grps[locus + group_name] = group_alleles_sort
                    if dev:
                        print('g-groups --> locus:', locus, ' -- group_alleles:', group_alleles_sort, ' -- group_name:', group_name)
    return g_grps


def build_dna_ser_table():
    dna_ser_tbl = {}
    with open(args.dstable) as f:
        for line in f:
            if not line.startswith("#"):
                line_strip = line.rstrip()
                line_values = line_strip.split(';')
                if line_values[2] != '' and line_values[2] != '?' and line_values[2] != '0':
                    locus = line_values[0].replace('*', '')  # locus without star symbol
                    allele = ":".join(line_values[1].split(":", 2)[:2])  # only 4 digit
                    dna_ser_tbl[locus + allele] = line_values[2]  # use locus + allele as dict key
                    if dev:
                        print('dna-ser -->  locus:', locus, ' -- allele:', allele, ' -- serology:', line_values[2])
    return dna_ser_tbl


def clean_g_group_alleles(g_group_alleles):
    group_alleles = set()
    alleles = g_group_alleles.split('/')
    for allele in alleles:
        if 'N' not in allele and 'Q' not in allele:
            clean_allele = ":".join(allele.split(":", 2)[:2])
            if clean_allele not in group_alleles:
                group_alleles.add(clean_allele)
    return group_alleles


def read_allelelist(allelelist, genotypes):
    frequent_alleles = {'A': set(), 'B': set(), 'C': set(), 'DRB1': set(), 'DQB1': set()}
    for genotype in genotypes['genotypes']:
        frequent_alleles["A"].update(genotype['A1'])
        frequent_alleles["A"].update(genotype['A2'])
        frequent_alleles["B"].update(genotype['B1'])
        frequent_alleles["B"].update(genotype['B2'])
        frequent_alleles["C"].update(genotype['C1'])
        frequent_alleles["C"].update(genotype['C2'])
        frequent_alleles["DRB1"].update(genotype['DRB11'])
        frequent_alleles["DRB1"].update(genotype['DRB12'])
        frequent_alleles["DQB1"].update(genotype['DQB11'])
        frequent_alleles["DQB1"].update(genotype['DQB12'])

    rare_alleles = {}
    with open(allelelist, "r") as input:
        for row in input:
            if not row.startswith("#"):
                cols = row.strip().split(",")
                allele = cols[1]
                if "*" in allele:
                    locus = allele.split("*")[0]
                    if locus in frequent_alleles and locus not in rare_alleles:
                        rare_alleles[locus] = set()
                    allele = allele.split("*")[1]
                    if ":" in allele:
                        allele = allele.split(":")[0] + ":" + allele.split(":")[1]
                        if locus in frequent_alleles and allele not in frequent_alleles[locus]:
                            rare_alleles[locus].add(allele)
    return rare_alleles

def parse_response_single(response_json, api_request_data, results):
    response_raw_p1 = response_json["pircheI"]
    response_raw_p2 = response_json["pircheII"]
    if "id_map" in api_request_data and len(response_raw_p1.keys()) > 0 and len(response_raw_p1.values()) and len(response_raw_p2.values()) > 0:
        if list(response_raw_p1.keys())[0] == api_request_data["id_map"]["d_id"][1]:
            result_data = {'id': api_request_data["id_map"]["d_id"][0], 'pircheI_scores': list(response_raw_p1.values())[0], 'pircheII_scores': list(response_raw_p2.values())[0]}
            if dev:
                print(result_data)
            results.append(result_data)
        else:
            print('ERROR: request(' + api_request_data["id_map"]["d_id"][1] + ') and response (' + list(response_raw_p1.keys())[0] + ') ids do not match. Result for donor_ID (' + api_request_data["id_map"]["d_id"][0] + ') skipped.')

def write_results_single(match_results):
    with open(args.output, "w") as output:
        writer = csv.writer(output, delimiter=',')
        writer.writerow(["donor_id", "PIRCHE_II", "PIRCHE_II_A", "PIRCHE_II_B", "PIRCHE_II_C", "PIRCHE_II_DRB1", "PIRCHE_II_DQB1"])

        for match_result in match_results:
            p2_scores = match_result['pircheII_scores']
            writer.writerow([str(match_result['id']), p2_scores['sum'], p2_scores['A'], p2_scores['B'], p2_scores['C'], p2_scores['DRB1'], p2_scores['DQB1']])


def parse_response_multi(response_text, id_list, results):
    data = {}
    mapping = {}
    header = 2
    try:
        lines = response_text.splitlines()
        for line in lines:
            cols = line.strip().split(",")
            if header == 1:
                for idx, k in enumerate(cols):
                    mapping[k] = idx
            if not header > 0:
                id = cols[0]
                if id in id_list:
                    real_don_id = id_list[id]["d_id"][0]
                    cols[0] = real_don_id
                    data[real_don_id] = cols
            header -= 1
    except UnicodeDecodeError:
        print("Skip")

    for key in data:

        new_row = []
        new_row.append(data[key][0])  # id
        new_row.append(data[key][28])  # pirche II sum

        new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_A_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_A_Epitopes"))
        new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_B_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_B_Epitopes"))
        new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_C_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_C_Epitopes"))
        new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DRB1_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DRB1_Epitopes"))
        new_row.append(get_column_score(data[key], mapping, "DRB1_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DRB3_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DRB4_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DRB5_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DQA1_DQB1_Presents_DQB1_Epitopes") + get_column_score(data[key], mapping, "DPA1_DPB1_Presents_DQB1_Epitopes"))

        results.append(new_row)


def get_column_score(cols, mapping, col_name):
    if col_name in mapping:
        score = 0
        col = cols[mapping[col_name]]
        peps = col.split("|")
        for pep in peps:
            if "[" in pep and not "{d}" in pep:
                s = pep.split("[")[1]
                s = float(s.split("]")[0])
                score += s
        return score
    return -1

def write_results_multi(results):

    with open(args.output, "w") as output:
        writer = csv.writer(output, delimiter=';')
        writer.writerow(["donor_id", "PIRCHE_II", "PIRCHE_II_A", "PIRCHE_II_B", "PIRCHE_II_C", "PIRCHE_II_DRB1", "PIRCHE_II_DQB1"])

        for result in results:
            writer.writerow(result)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-d", "--dev", help="dev mode", action="store_true")
    parser.add_argument("-url", "--url", help="target url", type=str, required=True)
    parser.add_argument("-u", "--user", help="username for request", type=str)
    parser.add_argument("-p", "--password", help="user password", type=str)
    parser.add_argument("-k", "--apikey", help="api key", type=str)
    parser.add_argument("-i", "--input", help="typing data input", type=str, required=True)
    parser.add_argument("-o", "--output", help="output file name", required=True)
    parser.add_argument("-a", "--anonymization", help="Enable anonymization. Default - True - anonymization enabled", default=True)
    parser.add_argument("-s", "--salt", help="Salt (password) used to anonymize input data. Use identical password when submitting same HLA data set multiple times.")
    parser.add_argument("-ka", "--kanonymization", help="Number of smoke hla data sets (genotypes) generated per patient and per donor. Default - 3.", default=5)
    parser.add_argument("-rm", "--requestmode", help="Single (single) or multiple (multi) patient donor pairs per request. Default - multi.", type=str, default="multi")
    parser.add_argument("-rs", "--requestsize", help="Max number of patient donor pairs per request for multi request setting. Max 500. Default - 300.", type=int, default=300)
    parser.add_argument("-pp", "--population", help="Population for HLA typing data provided (needed for low res high res conversion). Default - NMDP EUR haplotypes (2007).", type=str, default="NMDP EUR haplotypes (2007)")
    parser.add_argument("-gg", "--ggroups", help="HLA g-groups reference table file name and path. Default - hla_nom_g.txt ", default="hla_nom_g.txt")
    parser.add_argument("-ds", "--dstable", help="HLA dna ser reference table file and path. Default - rel_dna_ser.txt", default="rel_dna_ser.txt")
    parser.add_argument("-al", "--allelelist", help="HLA allelelist file. Default - Allelelist.txt", default="Allelelist.txt")
    parser.add_argument("-ht", "--haplotypes", help="NMDP haplotype table file (either 2007 or 2011 or equally formatted). Default - 2007_haplotypes.xls", default="2007_haplotypes.xls")
    parser.add_argument("-hf", "--haplofileformat", help="NMDP haplotype table file alleles format (either alleles XXXX (2007) or locus + alleles L*XX:XX (2011)). Default - XXXX", default="XXXX")
    parser.add_argument("-hp", "--haplofilepop", help="NMDP haplotype table file population short code as used in the header row (e.g. EUR [2007] or EURCAU [2011] ). Default - EUR", default="EUR")
    parser.add_argument("-hth", "--haplothreshold", help="frequency threshold for haplotypes 0.0 to 1.0. Default - 0.8", default="0.8")
    parser.add_argument("-prx", "--proxyhttp", help="HTTP proxy - full http proxy url (ip or dns) with protocol and port (http(s)://proxy:port)")
    parser.add_argument("-prxs", "--proxyhttps", help="HTTPS proxy - full https proxy url (ip or dns) with protocol and port (http(s)://proxy:port)")

    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        print("verbose mode active")

    dev = args.dev
    if dev:
        print("dev mode active")

    print("Initialize base data...")

    g_groups_global = build_g_groups()
    dna_ser_table_global = build_dna_ser_table()
    if args.anonymization:
        genotype_data = build_genotypes()
    else:
        genotype_data = []

    rare_alleles = read_allelelist(args.allelelist, genotype_data)

    raw_csv_data = {}

    referenceHeader = 'pid,pA1,pA2,pB1,pB2,pC1,pC2,pDRB11,pDRB12,pDRB31,pDRB32,pDRB41,pDRB42,pDRB51,pDRB52,' \
                      'pDQA11,pDQA12,pDQB11,pDQB12,pDPA11,pDPA12,pDPB11,pDPB12,' \
                      'did,dA1,dA2,dB1,dB2,dC1,dC2,dDRB11,dDRB12,dDRB31,dDRB32,dDRB41,dDRB42,dDRB51,dDRB52,' \
                      'dDQA11,dDQA12,dDQB11,dDQB12,dDPA11,dDPA12,dDPB11,dDPB12'

    print("Import and parse data...")
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

    proxies = {
        "http": args.proxyhttp,
        "https": args.proxyhttps
    }

    results = []

    if args.requestmode == "multi":
        print("Generating smoke data...")
        api_request_data = get_api_request_data_multi(raw_csv_data, genotype_data, proxies)
        print("Performing requests and parsing results...")
        numberRequests = str(len(api_request_data["api_payloads"]))
        requestCount = 1
        for api_payload in api_request_data["api_payloads"]:
            print("Running request #", requestCount, " of ", numberRequests)
            requestCount += 1
            response = rest_request(args.url, args.user, args.password, args.apikey, api_payload, proxies)

            if response.status_code == 200:
                parse_response_multi(response.text, api_request_data["id_list"], results)
            else:
                print('ERROR: a request failed')

        print("Writing results...")
        write_results_multi(results)

    elif args.requestmode == "single":
        print("Generating smoke data...")
        api_requests_data = get_api_requests_data_single(raw_csv_data, genotype_data)
        print("Performing requests and parsing results...")
        numberRequests = str(len(api_requests_data))
        requestCount = 1
        for api_request_data in api_requests_data:
            print("Running request #", requestCount, " of ", numberRequests)
            requestCount += 1
            response = rest_request(args.url, args.user, args.password, args.apikey, api_request_data["api_payload"], proxies)
            if response.status_code == 200:
                parse_response_single(response.json(), api_request_data, results)
            else:
                print('ERROR: a request failed')

        print("Writing results...")
        write_results_single(results)
