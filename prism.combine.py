#!/usr/bin/env python3

import re
import os
import argparse
import pandas as pd
import sqlite3
import warnings


def main():
    parser = argparse.ArgumentParser(description="PRISM combiner")

    parser.add_argument('-i', nargs='+', default=[], required=True, help="input CSV file list (files can be gzipped)")
    parser.add_argument('-s', required=True, help='tab delimited file with the description of samples: '
                                                  'Source_File	Sample_Name	Sample_Replica	Sample_Type')
    parser.add_argument('-p', nargs='+', default=[], required=False,
                        help='list of PRISM databases used for input files (prio1 or prio2 or extra or smth else)')
    parser.add_argument('-d', default='\t', required=False, help='input CSV file delimiter')
    parser.add_argument('-db', required=False, help='database file (if not specified, the data will store in memory)')
    parser.add_argument('-q', default=0.01, required=False, help='Q threshold (default Q<0.01 = 1%FDR)')
    parser.add_argument('-o', default='combined_results.csv.gz', required=False, help='output file with results')

    args = parser.parse_args()

    input_file_list = args.i
    sample_file = args.s
    output_file = args.o
    run_names = args.p
    csv_delimiter = args.d
    db_file = args.db
    q_threshold = args.q

    if not len(run_names):
        for file_name in input_file_list:
            file_name = re.sub(r'.*\/', '', file_name)
            file_name = re.sub(r'\..*', '', file_name)
            run_names.append(file_name)
    else:
        if len(run_names) != len(input_file_list):
            print("Wrong amount of run names: {} (should be {})".format(len(run_names), len(input_file_list)))
            exit(1)

    all_data = []
    for i, file_name in enumerate(input_file_list):
        if re.search(r'\.gz$', file_name):
            data = pd.read_csv(file_name, compression='gzip', sep=csv_delimiter, quotechar='"', low_memory=False)
        else:
            data = pd.read_csv(file_name, sep=csv_delimiter, quotechar='"', low_memory=False)

        data['Databases_PRISM'] = run_names[i]
        all_data.append(data)

    all_data = pd.concat(all_data)

    warnings.simplefilter(action='ignore', category=FutureWarning)  # to suppress FutureWarning
    all_data.columns = all_data.columns.str.replace('[%()/*:]', '')
    all_data.columns = all_data.columns.str.strip().str.replace('[ .-]', '_')
    all_data.columns = all_data.columns.str.strip().str.replace('_+', '_')

    if db_file:
        if os.path.exists(db_file):
            os.remove(db_file)
        connector = sqlite3.connect(db_file)
    else:
        connector = sqlite3.connect(':memory:')

    sample_description = pd.read_csv(sample_file, sep="\t", quotechar='"', low_memory=False)
    sample_description.to_sql(name='description', con=connector, index=False)

    all_data.index += 1  # to start the index with 1, but not zero
    all_data.index.name = 'Rec_ID'
    all_data.to_sql(name='data', con=connector)

    cur = connector.cursor()
    sql = 'CREATE INDEX source_file_index ON description(Source_File);'
    cur.execute(sql)

    sql = 'CREATE TABLE peptides AS SELECT DISTINCT Sequence FROM data;'
    cur.execute(sql)
    connector.commit()

    sql = 'CREATE TABLE ext_data AS SELECT * FROM data INNER JOIN description USING(Source_File);'
    cur.execute(sql)
    connector.commit()

    sql = 'CREATE INDEX sequence_index ON ext_data(Sequence);'
    cur.execute(sql)

    sql = 'SELECT COUNT(*) AS n FROM ext_data;'
    n_rec_ext_data = pd.read_sql(sql, connector).iloc[0]['n']
    sql = 'SELECT COUNT(*) AS n FROM data;'
    n_rec_data = pd.read_sql(sql, connector).iloc[0]['n']

    if n_rec_ext_data != n_rec_data:
        print("You probably have an incorrect file with the description of the samples")
        exit(1)

    # FROM Andreas's email:
    # I typically filter by Q (e.g. Q<0.01 = 1%FDR) and by NetMHC predictions (e.g. rank<2%).
    # I generate this column in Spotfire. Best Q is the lowest Q value for a peptide sequence among all merged samples.
    # I.e., when a peptide is reliably identified in sample 1 (e.g. Q<0.01) and in sample 2 Q>0.01,
    # the peptide is still counted as identified in both samples. This increases the overlap between merged samples.

    data_output = None

    sql = 'SELECT Sequence FROM peptides;'
    peptides = pd.read_sql(sql, connector)
    for peptide in peptides['Sequence']:
        print(peptide)
        sql = 'SELECT * FROM ext_data WHERE Sequence = "{}" ORDER BY Q, netMHC_rank LIMIT 1;'.format(peptide)
        best_peptide_rec = pd.read_sql(sql, connector)
        if len(best_peptide_rec.index) and best_peptide_rec.iloc[0]['Q'] < q_threshold:
            best_q = best_peptide_rec.iloc[0]['Q']  # minimum false discovery rate (FDR)

            sql = 'SELECT ALC FROM ext_data WHERE Sequence = "{}" ORDER BY ALC DESC LIMIT 1;'.format(peptide)
            best_alc_rec = pd.read_sql(sql, connector)
            best_alc = best_alc_rec.iloc[0]['ALC']  # average local confidence (ALC)

            best_peptide_rec['HLA_allele'] = best_peptide_rec['HLA_allele'].str.replace(r'[*:]', '')
            best_peptide_rec['HLA_allele'] = best_peptide_rec['HLA_allele'].str.replace(r'-', '_')

            best_hla = best_peptide_rec.iloc[0]['HLA_allele']
            filtered_hla = best_hla if best_hla in best_peptide_rec.columns and best_peptide_rec.iloc[0][best_hla] < 2 else ''

            sql = 'SELECT group_concat(DISTINCT Category) AS Categories FROM ext_data WHERE Sequence = "{}";'.format(peptide)
            categories_rec = pd.read_sql(sql, connector)
            categories = categories_rec.iloc[0]['Categories']

            sql = 'SELECT group_concat(DISTINCT Sample_Type) AS Types FROM ext_data ' \
                  'WHERE Sequence = "{}";'.format(peptide)
            status_rec = pd.read_sql(sql, connector)
            status = status_rec.iloc[0]['Types']

            sql = 'SELECT group_concat(DISTINCT Databases_PRISM) AS Databases_PRISMs FROM ext_data ' \
                  'WHERE Sequence = "{}";'.format(peptide)
            db_prism_rec = pd.read_sql(sql, connector)
            db_prism = db_prism_rec.iloc[0]['Databases_PRISMs']

            sql = 'SELECT DISTINCT Sample_Name FROM ext_data WHERE Sequence = "{}";'.format(peptide)
            sample_rec = pd.read_sql(sql, connector)
            sample_array = sample_rec['Sample_Name'].values
            samples = ','.join(sample_array)

            sql = 'SELECT DISTINCT Sample_Name FROM description ORDER BY Sample_Name;'
            samples_rec = pd.read_sql(sql, connector)
            intensity_sum = 0
            intensities = []
            replica_counts = []
            for sam_index, sam_row in samples_rec.iterrows():
                sample_name = sam_row['Sample_Name']
                sql = 'SELECT DISTINCT Sample_Replica FROM ext_data ' \
                      'WHERE Sequence = "{}" AND Sample_Name = "{}";'.format(peptide, sample_name)
                replica_count_rec = pd.read_sql(sql, connector)
                replica_counts.append({'sample': sample_name, 'rep_count': len(replica_count_rec['Sample_Replica'].values)})

                sql = 'SELECT DISTINCT Sample_Replica FROM description ' \
                      'WHERE Sample_Name = "{}" ORDER BY Sample_Replica;'.format(sample_name)
                replica_rec = pd.read_sql(sql, connector)
                for rep_index, rep_row in replica_rec.iterrows():
                    sample_replica = rep_row['Sample_Replica']
                    sql = 'SELECT Intensity FROM ext_data ' \
                          'WHERE Sequence = "{}" AND Sample_Name = "{}" AND Sample_Replica = "{}" LIMIT 1' \
                          ';'.format(peptide, sample_name, sample_replica)
                    intensity_rec = pd.read_sql(sql, connector)
                    intensity = intensity_rec.iloc[0]['Intensity'] if len(intensity_rec.index) else 0

                    intensities.append({'replica': '_'.join([sample_name, sample_replica]), 'intensity': intensity})
                    intensity_sum += intensity

            hla_columns = [col for col in best_peptide_rec.columns if re.match(r'HLA_[ABC]', col)]
            base_columns = ['Source_File', 'Feature', 'Scan', 'ALC', 'Length', 'RT', 'Mass', 'ppm', 'ID',
                            'Location_count', 'Genome', 'Location', 'Sequence', 'Top_location_count',
                            'Top_location_count_no_decoy', 'Q', 'Gene', 'Symbol', 'ORF_location', 'HLA_allele',
                            'netMHC_rank']
            data_output_rec = best_peptide_rec[base_columns + hla_columns]
            data_output_rec['Filtered_HLA_allele'] = filtered_hla
            data_output_rec['Best_Q'] = best_q
            data_output_rec['Best_ALC'] = best_alc
            data_output_rec['Categories'] = categories
            data_output_rec['Status_over_sequence'] = status
            data_output_rec['Databases_PRISM'] = db_prism
            data_output_rec['Samples'] = samples
            for index, item in enumerate(replica_counts):
                data_output_rec[item['sample']] = item['rep_count']
            data_output_rec['Intensity_Sum'] = intensity_sum
            for index, item in enumerate(intensities):
                data_output_rec['Intensity_' + intensities[index]['replica']] = intensities[index]['intensity']
            # data_output_rec['Disease_states_over_ORF'] = '?'  # TODO ask Andreas about Disease_states_over_ORF

            if type(data_output) is pd.DataFrame:
                data_output = pd.concat([data_output, data_output_rec])
            else:
                data_output = data_output_rec

    if re.search(r'\.gz$', output_file):
        data_output.to_csv(output_file, compression='gzip', sep='\t', index=False)
    else:
        data_output.to_csv(output_file, sep='\t', index=False)

    print("...done")

# end of main()


if __name__ == '__main__':
    main()
