import gffutils

def extract_introns(gff_file, gff_db, bed_file):
    print("gff_file: ", gff_file)
    print("gff_db: ", gff_db)
    print("bed_file: ", bed_file)
    
    # Create the GFF database
    db = gffutils.create_db(gff_file, dbfn=gff_db, force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True, verbose=True)

    # Load the GFF database
    # db = gffutils.FeatureDB('gff_db', keep_order=True)

    genes = db.features_of_type("gene")

    junc_count = 0 
    with open(bed_file, 'w') as fw:
        for gene in genes:
            # Retrieve child features
            children = db.children(gene, level=1)

            # Print child features
            for child in children:
                exons = db.children(child, featuretype='exon', order_by='start', level=1)
                prev_exon = ""
                for exon in exons:

                    if (prev_exon != ""):
                        # print(exon[0] + "\t" + str(prev_exon.end-1) + "\t" + str(exon.start))
                        fw.write(exon[0] + "\t" + str(prev_exon.end-1) + "\t" + str(exon.start) + "\t" + f'JUNC{junc_count:08}' + "\t.\t" + exon.strand +"\n")
                        junc_count += 1

                    prev_exon = exon