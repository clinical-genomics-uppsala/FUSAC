def mol_count(input_dict):
    ref_sup = 0
    mut_sup = 0
    ffpe_sup = 0
    n_sup = 0
    del_sup = 0
    # Counts the total support for each variant-option for every UMI-key entry in input_dict
    for umi_key in input_dict:
        if input_dict[umi_key]["Mate_Hits"]:
            ref_sup += input_dict[umi_key]["Mate_Hits"]["Reference_Support"]
            mut_sup += input_dict[umi_key]["Mate_Hits"]["Mutation_Support"]
            ffpe_sup += input_dict[umi_key]["Mate_Hits"]["FFPE_Support"]
            n_sup += input_dict[umi_key]["Mate_Hits"]["N_Support"]
            del_sup += input_dict[umi_key]["Mate_Hits"]["Del_Support"]
        else:
            continue
    return [ref_sup, mut_sup, ffpe_sup, n_sup, del_sup]


def sup_count(input_dict, nuc):
    n_sup = {"Paired": {}, "String_1_Single": {}, "String_2_Single": {}}
    for n in nuc:
        n_sup["Paired"] = {"String_1": {}, "String_2": {}}
        n_sup["Paired"]["String_1"][n] = 0
        n_sup["Paired"]["String_2"][n] = 0
        n_sup["String_1_Single"][n] = 0
        n_sup["String_2_Single"][n] = 0
        for umi_key in input_dict:
            # Iterates through each alternative allele called by the variant caller
            # Counts the no. molecules supporting the variant for each direction in the single hits
            if input_dict[umi_key]:
                if input_dict[umi_key]["Single_Hits"]:
                    for sh_dir in input_dict[umi_key]["Single_Hits"]:
                        if input_dict[umi_key]["Single_Hits"][sh_dir]:
                            if n in input_dict[umi_key]["Single_Hits"][sh_dir]:
                                n_sup[sh_dir][n] += 1
            # Sees if the dict "Mate Hits" is empty, if not, counts the no molecules supporting the variant for each
            # dict within the "Mate Hits" dictionary
            if input_dict[umi_key]["Mate_Hits"]:
                new_lst = list(input_dict[umi_key]["Mate_Hits"].items())
                new_dict = dict(new_lst[0:-7])
                for v_h in new_dict:
                    if new_dict[v_h]:
                        for s_t in new_dict[v_h]:
                            if new_dict[v_h][s_t]:
                                if n in new_dict[v_h][s_t]:
                                    n_sup["Paired"][s_t][n] += 1
    return n_sup
