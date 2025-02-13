#!/usr/bin/python3

import json

dir_proj = "/mnt/chaelab/rachelle/zzOtherzz/XiaoMei/HDstress"

networks = {"TFTF": "TFTF_simplify", "keyGRN": "keyGRN_simplify", "GRN": "GRN_simplify"}
network_name = "TFTF"

for network_name in networks:
    ## parse cx2 file
    f_cx2 = dir_proj + f"/data/network/{networks[network_name]}.edges.cx2"
    with open(f_cx2, 'r') as f:
        data = json.load(f)
        dat_nodes = [x for x in data if "nodes" in x]
        dat_nodes = dat_nodes[0]["nodes"]
    ## write
    fout = f_cx2 + ".cytoscape-pfd.coord"
    with open(fout, "w+") as f:
        _ = f.write('\t'.join(["node", "x", "y"]) + '\n')
        for node in dat_nodes:
            _ = f.write('\t'.join(map(str, [node["v"]["name"], node["x"], node["y"]])) + '\n')

