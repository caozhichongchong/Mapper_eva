import os,glob
import argparse
import pandas as pd

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="input sam file",
                      type=str, default='input.sam',
                      metavar='input.sam')

################################################## Definition ########################################################
# setting match score as 0 for all
args = parser.parse_args()
# AS_col = {
#     'bowtie':11,
#     'bwa':13,
#     'minimap':13,
#     'mapper': 11
# }
def samfilter(samfile):
    # for tool in AS_col:
    #     if tool in os.path.basename(samfile):
    #         break
    # AS_col_sam = AS_col[tool]
    for lines in open(samfile,'r'):
        if not lines.startswith('@') and not lines.startswith('#'):
            lines_set = lines.split('\t')
            for element in lines_set:
                if element.startswith('cs'): #mapper only
                    AS_col_sam=lines_set.index(element)
                    break
                elif element.startswith('AS'):
                    AS_col_sam=lines_set.index(element)
                    break
            break
    print(samfile,AS_col_sam,element)
    allsam = pd.read_csv(samfile,
                         sep='\t', usecols=[0, 1, 2, 3,AS_col_sam], header=None,
                         names=['readID', 'Direction', 'CHR', 'POS','AS'], comment='@')
    allsam = allsam[~allsam['AS'].isna()]
    allsam['score']=[float(x.split(':')[-1]) for x in allsam['AS']]
    # allsam = allsam.sort_values(by=['readID', 'score'], ascending=[True, False])
    # allsam['readIDAS']=allsam['readID'] + [str(x) for x in allsam['AS']]
    # allsamkeep = allsam.drop_duplicates(['readID','Direction'])
    # allsam = allsam[allsam['readIDAS'].isin(allsamkeep['readIDAS'])]
    allsam = allsam.drop_duplicates(['readID', 'Direction','CHR', 'POS','AS','score'])
    allsam.loc[:,['readID','Direction', 'CHR', 'POS','AS','score']].to_csv(samfile,sep='\t',index=False)

samfilter(args.i)

