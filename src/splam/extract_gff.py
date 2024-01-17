import json
import gffutils

def write_introns(fw, intron_dict, junc_count):
    for id, intron_info in intron_dict.items():
        eles = id.split("\t")
        eles.insert(3, intron_info["count"])
        eles.append(intron_info["trans"])
        # print("\t".join(eles))
        trans_ls = ",".join(eles[5])
        fw.write(f'{eles[0]}\t{str(eles[1])}\t{str(eles[2])}\tJUNC{junc_count:08}\t{str(eles[3])}\t{str(eles[4])}\t{trans_ls}\n') 
        junc_count += 1
    intron_dict = {}
    return junc_count, intron_dict

def get_parent_features_to_extract(feature_types_file):
    feature_types = ["gene"]
    if feature_types_file is not None:
        feature_types = []
        f = open(feature_types_file)
        for line in f.readlines():
            feature_types.append(line.rstrip())
    return feature_types


def extract_introns(gff_file, gff_db, feature_file, is_load_gff_db, bed_file, trans_intron_num_txt):
    print("gff_file: ", gff_file)
    print("gff_db: ", gff_db)
    print("bed_file: ", bed_file)
    
    if not is_load_gff_db:
        # Create the GFF database
        db = gffutils.create_db(gff_file, dbfn=gff_db, force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True, verbose=True)
    else:
        # Load the GFF database
        db = gffutils.FeatureDB(gff_db, keep_order=True)

    features = get_parent_features_to_extract(feature_file)
    intron_dict = {}
    trans_2_intron_num = {}
    junc_count = 0 
    with open(bed_file, 'w') as fw:
        bundle_s = 0
        bundle_e = 0
        bundle_chr = ""
        for feature in features:
            for locus in db.features_of_type(feature):
                bundle_s, bundle_e, bundle_chr, junc_count, intron_dict= extract_introns_itr(db, locus, bundle_s, bundle_e, bundle_chr, intron_dict, trans_2_intron_num, junc_count, fw)

        # Write out all introns
        junc_count, intron_dict = write_introns(fw, intron_dict, junc_count)
        intron_dict = {}

    json.dump(trans_2_intron_num, open(trans_intron_num_txt,'w'))



def extract_introns_itr(db, locus, bundle_s, bundle_e, bundle_chr, intron_dict, trans_2_intron_num, junc_count, fw):
    children = list(db.children(locus, featuretype='exon', level=1, order_by='start'))
    if len(children) == 0:
        children = db.children(locus, level=1)
        for child in list(children):
            bundle_s, bundle_e, bundle_chr, junc_count, intron_dict = extract_introns_itr(db, child, bundle_s, bundle_e, bundle_chr, intron_dict, trans_2_intron_num, junc_count, fw)

    else:
        if bundle_chr != locus[0]:
            # Reconstruct the bundle
            bundle_s = locus.start
            bundle_e = locus.end
            bundle_chr = locus[0]
            # Write out all introns
            junc_count, intron_dict = write_introns(fw, intron_dict, junc_count)

        elif bundle_chr == locus[0] and locus.start <= bundle_e:
            # Extend the bundle
            bundle_e = locus.end

        elif bundle_chr == locus[0] and locus.start > bundle_e:
            bundle_s = locus.start
            bundle_e = locus.end
            # Write out all introns
            junc_count, intron_dict = write_introns(fw, intron_dict, junc_count)

        # Print child features
        exons = db.children(locus, featuretype='exon', level=1, order_by='start')
        prev_exon = ""
        for exon in exons:
            if (prev_exon != ""):
                key = f'{exon[0]}\t{str(prev_exon.end-1)}\t{str(exon.start)}\t{exon.strand}'
                if locus.id not in trans_2_intron_num.keys():
                    trans_2_intron_num[locus.id] = 1    
                else:
                    trans_2_intron_num[locus.id] += 1
                if key in intron_dict.keys():
                    intron_dict[key]["count"] += 1
                    intron_dict[key]["trans"].add(locus.id)
                else:
                    intron_dict[key] = {}
                    intron_dict[key]["count"] = 1
                    intron_dict[key]["trans"] = set([locus.id])
            prev_exon = exon
    return bundle_s, bundle_e, bundle_chr, junc_count, intron_dict

    # with open(bed_file, 'w') as fw:

    #     bundle_s = 0
    #     bundle_e = 0
    #     bundle_chr = ""

    #     for gene in genes:

    #         # curr_s = gene.start
    #         # curr_e = gene.end

    #         if bundle_chr != gene[0]:
    #             # Reconstruct the bundle
    #             bundle_s = gene.start
    #             bundle_e = gene.end
    #             bundle_chr = gene[0]
    #             # Write out all introns
    #             junc_count = write_introns(fw, intron_dict, junc_count)
    #             intron_dict = {}

    #         elif bundle_chr == gene[0] and gene.start < bundle_e:
    #             # Extend the bundle
    #             bundle_e = gene.end

    #         # Retrieve child features
    #         children = db.children(gene, level=1)

    #         # Print child features
    #         for child in children:
    #             chrom = child[0]
    #             trans_id = child.id
    #             trans_s = child.start
    #             trans_e = child.end
    #             trans_2_intron_num[trans_id] = 0

    #             exons = db.children(child, featuretype='exon', order_by='start', level=1)
    #             prev_exon = ""
    #             for exon in exons:
    #                 if (prev_exon != ""):
    #                     key = f'{exon[0]}\t{str(prev_exon.end-1)}\t{str(exon.start)}\t{exon.strand}'
    #                     trans_2_intron_num[trans_id] += 1
    #                     if key in intron_dict.keys():
    #                         intron_dict[key]["count"] += 1
    #                         intron_dict[key]["trans"].add(trans_id)
    #                     else:
    #                         intron_dict[key] = {}
    #                         intron_dict[key]["count"] = 1
    #                         intron_dict[key]["trans"] = set([trans_id])
    #                 prev_exon = exon
        
    #     # Write out all introns
    #     junc_count = write_introns(fw, intron_dict, junc_count)
    #     intron_dict = {}
    
    # json.dump(trans_2_intron_num, open(trans_intron_num_txt,'w'))