
import argparse
import json
import pandas as pd


def group_lineages(tsv_file, lin_list, other_lins, alias, lin_map, other_map, threshold=0.01):
    results = {}

    lin_file = pd.read_csv(tsv_file, sep='\t', index_col=0).to_dict()
    values = lin_file[list(lin_file.keys())[0]]

    # Empty results file
    if values['summarized'] == '[]':
        return

    for lineage, abundance in zip(values['lineages'].split(), values['abundances'].split()):
        if lineage[0] != 'X':
            prefix = lineage.split('.')
            prefix[0] = alias[prefix[0]]
            expanded = '.'.join(prefix)
            # Check if the expanded form is a particular lineage being looked for
            if expanded in lin_list:
                group = lin_map[expanded]
            else:
                # Try grouping it to the nearest ancestor
                # other_lins is a reverse sorted list. e.g. ['recombinant', 'B.1.1.529.5.3.1.1.1.1', 'B.1.1.529.5.2.1', 'B.1.1.529.5']
                for lin in other_lins:
                    if lin == 'recombinant':
                        continue
                    if lin in expanded:
                        group = other_map[lin]
                        break
                if group == None:
                    group = 'Other'
        else:
            # Recombinant
            if lineage in lin_list:
                group = lin_map[lineage]
            else:
                group = 'Other recombinants'
                for lin in other_lins:
                    if lin in lineage:
                        group = other_map[lin]
                        break

        results.update({lineage: {'abundance': float(abundance), 'ancestor': group if float(abundance) > threshold else "Minor"}})
    
    return results


def add_ci(sample_result, boot_file):
    result_file = pd.read_csv(boot_file, sep=',', index_col=0).to_dict()
    for lineage in sample_result.keys():
        sample_result[lineage].update({'lower.95' : result_file[lineage][0.025], 'upper.95' : result_file[lineage][0.975]})

    return sample_result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate summary file"
    )
    parser.add_argument('linfile', type=str,
                        help="Path to lin TSV")
    parser.add_argument('bootfile', type=str,
                        help="Path to the bootstrap lineage CSV")   
    parser.add_argument('lab', type=str,
                        help="Which lab (western/guelph/waterloo)")
    parser.add_argument('run', type=str,
                        help="Run Name")   
    parser.add_argument('sample', type=str,
                        help="Sample Name")         
    parser.add_argument('--alias', type=str, default="data/alias_key.json",
                        help="<input> PANGO aliases")
    args = parser.parse_args()

    with open(args.alias, 'r') as alias_file:
        alias = json.loads(alias_file.read())

    # Lineages of interest
    loi = ['Other recombinants', 'XBB.1.5', 'XBB.1', 'XBM', 'Other BA.5', 'Other BA.4', 'Other BA.2', 'Other BA.1', 'Other BQ', 'Other BF', 'BQ.1.1', 'BQ.1', 'BQ.1.5', 'BF.7', 'BQ.1.1.1', 'BQ.1.1.4', 'CH.1.1', 'BQ.1.13', 'BA.5.2', 'BA.5.2.1']
    
    # Map lineages to expanded form if they are not a recombinant e.g. B.1.1.529.5.3.1.1.1.1.1 : BQ.1
    lin_map = {}

    # Map other lineages to expanded form if they are not a recombinant e.g. B.1.1.529.5 : "Other BA.5"
    other_map = {}

    for lin in loi:

        # Split by space in case lineage of interest is "Other BA.5" for example
        lin_string = lin.split() 
        
        # Handle Case when there are "Other" lineages of interest
        if lin_string[0] == 'Other':
            prefix = lin_string[1].split('.')
            if lin_string[1][0] == 'X':
                # Case when we want to group "Other XBB.1" for example
                other_map.update({lin_string[1] : lin_string[1] })
                continue

            try:
                prefix[0] = alias[prefix[0]]
                expanded = '.'.join(prefix)
                other_map.update({expanded : lin})
            except KeyError:
                # KeyError for "Other recombinants"
                other_map.update({'recombinant' : lin})

        # If the lineage of interest is a recombinant, do not look into its alias
        elif lin[0] == 'X':
            lin_map.update({lin : lin})
        else:
            prefix = lin.split('.')
            prefix[0] = alias[prefix[0]]
            expanded = '.'.join(prefix)
            lin_map.update({expanded : lin})

    lin_list = sorted(lin_map.keys(), reverse=True)    
    # other_lins = sorted(other_map.keys(), reverse=True)
    other_lins = other_map.keys()

    sample_result = group_lineages(args.linfile, lin_list, other_lins, alias, lin_map, other_map)
    sample_result = add_ci(sample_result, args.bootfile)

    lineages = sorted(sample_result.keys())
    columns = ['estimated frequency', 'lower.95', 'upper.95', 'nearest ancestor (of interest)']

    df = pd.DataFrame(columns=columns, index=lineages)
    for lineage, values in sample_result.items():
        df['estimated frequency'][lineage] = sample_result[lineage]['abundance']
        df['lower.95'][lineage] = sample_result[lineage]['lower.95']
        df['upper.95'][lineage] = sample_result[lineage]['upper.95']
        df['nearest ancestor (of interest)'][lineage] = sample_result[lineage]['ancestor']

    df.to_csv('{}/{}/{}.freyja.csv'.format(args.lab, args.run, args.sample), encoding='utf-8')