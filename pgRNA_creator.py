#!/usr/bin/env python3

import re
import csv
import json
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats, integrate
import seaborn as sns


# modify annotation and include 5kb window to UTR
# the UTR annotations are added to the end of the GTF file
# if gene on plus strand then subtract 5000bp from Start site
# if gene on minus strand then add 5000bp to the Start site
def modify_anno(raw, mod):

    # Read annotation file and identify UTR region
    with open(raw, 'r') as r:

        # Write modified UTR at end of file
        with open(mod, 'a') as a:

            # parse raw line by line
            for line in r:

                split = line.split('\t')

                # check if feature is for UTR and on positive strand
                if split[7] == 'UTR' and split[5] == '+':

                    split[1] = str(int(split[1]) - 5000)

                    # ensure UTR doesn't fall off genome
                    if int(split[1]) < 0:

                        split[1] = '0'

                    # append updated annotation to file
                    a.write('\t'.join(split))

                # check if feature is for UTR and on negative strand
                elif split[7] == 'UTR' and split[5] == '-':

                    split[2] = str(int(split[2]) + 5000)

                    # append updated annotation to file
                    a.write('\t'.join(split))


def pgRNA_creator(raw, filter1):

    # Read 'GuideScan_batch_output.csv' and obtain workable sgRNAs
    with open(raw, 'r') as r:

        # file to store filtered sgRNAs
        with open(filter1, 'w') as w, open('onegRNA.txt', 'w') as s, open('pgRNA.txt', 'w') as p:
            
            # set target
            target = None

            # start gRNA count upon target match
            lcount = 0
            rcount = 0
            
            # parse raw line by line
            for line in r:
                
                # if line contains chr#:#-# then handle
                if re.match('chr\S+:\d+-\d+', line):
        
                    # check if target is None
                    if target is None:
                        
                        # target None so save information from the line
                        target = line
                        left = dict()
                        right = dict()

                    # if one pgRNA then write to onepgRNA.txt
                    elif (lcount and rcount) !=0 and (lcount + rcount) < 3:

                        # write information from previous target
                        s.write(target + '\n')
                        s.write('LEFT\n')
                        s.write(json.dumps(left) + '\n')
                        s.write('RIGHT\n')
                        s.write(json.dumps(right) + '\n\n')
                        
                    # otherwise write target, left & right information to file
                    elif (lcount and rcount) != 0 and (lcount + rcount) > 3:

                        # use Cartesian product for pgRNA Creation
                        product = list(itertools.product(left, right))

                        p.write('\n' + target + '\n')

                        v = 0
                        
                        # iterate and identify pgRNAs
                        while v < len(product):

                            # add cutting specificity score
                            sgRNA_left = product[v][0]
                            sgRNA_right = product[v][1]
                            seql = left[sgRNA_left][3]
                            seqr = right[sgRNA_right][3]

                            # ensure specificity values are present
                            try:

                                scoreCom = float(left[sgRNA_left][5]) + float(right[sgRNA_right][5])
                                specifCom = float(left[sgRNA_left][4]) + float(right[sgRNA_right][4])
                                totalCoord = int(right[sgRNA_right][2]) - int(left[sgRNA_left][1])

                                # write pgRNA to file
                                pgRNA_list = [sgRNA_left, seql, sgRNA_right, seqr, scoreCom, specifCom, abs(totalCoord)]
                                p.write(json.dumps(pgRNA_list) + '\n')
                                v += 1

                            except ValueError:

                                v += 1

                                pass
                        
                        # write information from previous target
                        w.write(target + '\n')
                        w.write('LEFT\n')
                        w.write(json.dumps(left) + '\n')
                        w.write('RIGHT\n')
                        w.write(json.dumps(right) + '\n\n')
                        
                        # store new target information
                        target = line
                        left = dict()
                        right = dict()

                        lcount = 0
                        rcount = 0
                        
                        # skip header line
                        next(r)

                    else:

                        target = line
                        left = dict()
                        right = dict()

                        lcount = 0
                        rcount = 0

                        next(r)
                        
                # if line contains chr#,# then handle
                if re.match('chr\S+,\d+', line):
                    
                    # split line by ',' and create list
                    sgRNA = line.split(',')

                    # handle index out of range error
                    try:
                        
                        # if column 10 contains annotation
                        if sgRNA[10] == '*' and lcount < 3:
                            
                            # save contents to left dictionary
                            left[sgRNA[11]] = sgRNA[0:8]
                            
                            # increment left count by 1
                            lcount += 1
                            
                    except IndexError:

                        pass

                    # handle index out of range error
                    try:
                        
                        # if column 24 contains annotation
                        if sgRNA[24] == '*' and rcount < 3:
                            
                            # save contents to right dictionary
                            right[sgRNA[25]] = sgRNA[14:22]
                            
                            # increment right count by 1
                            rcount += 1

                    except IndexError:

                        pass

def pgRNA_clean(raw, filter1):

    # final total gRNA and target sites
    ftotal = 0
    ftarget = 0

    # Open csv using CSV module
    with open(raw, newline='') as r:

        # file to store filtered sgRNAs
        with open(filter1, 'w') as f1, open('pgRNA.txt', 'w') as p:

            # set target
            target = None

            # start gRNA count upon target match
            lcount = 0
            rcount = 0

            # begin parsing CSV file
            reader = csv.reader(r)

            # row by row
            for row in reader:

                # if line contains chr#:#-# then handle 
                if re.match('chr\S+:\d+-\d+', ' '.join(row)):

                    # None means the first chr target is found
                    if target is None:

                        # target None so save information from the line
                        target = ' '.join(row)
                        left = dict()
                        right = dict()

                    # if left and right have gRNAs
                    # write information to file
                    elif (lcount and rcount) !=0:

                        # re-order according to specificity score
                        sleft = sorted(left.items(), key=lambda e: e[1][5], reverse = True)
                        sright = sorted(right.items(), key=lambda e: e[1][5], reverse = True)

                        # select top3 lgRNAs and rgRNAs w/ specificity > 0.7 & activity > 50
                        fleft, fright = finalClean(sleft, sright)

                        if fleft == False:
                            
                            continue

                        else:

                            # Set variable to test is target passes
                            target_pass = False

                            # use Cartesian product for pgRNA Creation
                            product = list(itertools.product(fleft, fright))

                            v = 0

                            # iterate and identify pgRNAs
                            while v < len(product):

                                # add cutting specificty score
                                sgRNA_left = product[v][0]
                                sgRNA_right = product[v][1]
                                seql = fleft[sgRNA_left][3]
                                seqr = fright[sgRNA_right][3]

                                # ensure specificity values are present
                                try:

                                    specificityCom = float(fleft[sgRNA_left][5]) + float(fright[sgRNA_right][5])
                                    activityCom = int(fleft[sgRNA_left][4]) + int(fright[sgRNA_right][4])
                                    coordCom = abs(int(fright[sgRNA_right][2]) - int(fleft[sgRNA_left][1]))

                                    if specificityCom >= 0.8 and target_pass == False:

                                        ftarget += 1
                                        target_pass = True
                                        p.write('\n' + target + '\n\n')
                                        
                                        # write pgRNA to file
                                        pgRNA_list = [sgRNA_left, seql, sgRNA_right, seqr, specificityCom, activityCom,
                                                      coordCom]
                                        p.write(json.dumps(pgRNA_list) + '\n')
                                        v += 1
                                        ftotal += 1

                                    elif specificityCom >= 0.8 and target_pass == True:

                                        # write pgRNA to file
                                        pgRNA_list = [sgRNA_left, seql, sgRNA_right, seqr, specificityCom, activityCom,
                                                      coordCom]
                                        p.write(json.dumps(pgRNA_list) + '\n')
                                        v += 1
                                        ftotal += 1

                                    else:

                                        v += 1
                                        pass

                                except ValueError:

                                    v += 1

                                    pass
                            
                            # write information from previous target
                            f1.write(target + '\n')
                            f1.write('LEFT\n')
                            f1.write(json.dumps(fleft) + '\n')
                            f1.write('RIGHT\n')
                            f1.write(json.dumps(fright) + '\n\n')

                        # store new target information
                        target = ' '.join(row)
                        left = dict()
                        right = dict()

                        lcount = 0
                        rcount = 0

                        next(r)

                # if line contains chr#,# then handle
                if re.match('chr\S+,\d+', ','.join(row)):

                    # handle index out of range error
                    try:

                        # save right gRNAs if not in coding region
                        if row[10] == '*':

                            # save contents to left dictionary
                            left[row[11]] = row[0:8]

                            # increment left count by 1
                            lcount += 1

                    except IndexError:

                        pass

                    # handle index out of range error
                    try:

                        # save left gRNAs if not in coding region
                        if row[24] == '*':

                            # save contents to right dictionary
                            right[row[25]] = row[14:22]

                            # increment right count by 1
                            rcount += 1

                    except IndexError:

                        pass

    print('Total Targets w/ gRNAs: ' + str(ftarget))
    print('Total gRNAs: ' + str(ftotal))
                        
# function to select top3 guides
# assuming specificity > 0.7 & activity > 50
def finalClean(sleft, sright):

    fleft = dict()
    fright = dict()
    tguide = 0
    lguide = 0
    rguide = 0

    # parse dictionary to identify top 3 left and right gRNAs
    for k in sleft:

        # ensure values are present
        try:
            if float(k[1][5]) >= 0.4 and int(k[1][4]) >= 50 and lguide < 3:
                
                fleft[k[0]] = k[1]
                lguide += 1
                tguide += 1
                
            elif lguide > 2:
                break

        except ValueError:

            if lguide > 2:
                break

    for k in sright:

        # ensure values are present
        try:
            
            if float(k[1][5]) >= 0.4 and int(k[1][4]) >= 50 and tguide < 6:
                
                fright[k[0]] = k[1]
                rguide += 1
                tguide += 1
                
            elif tguide > 5:
                break

        except ValueError:

            if tguide > 5:
                break

    if lguide and rguide > 0:

        return fleft, fright

    else:

        fleft, fright = False, False
        
        return fleft, fright


# parse pgRNA.txt and average coordiante distance
# for each of the pgRNA targets
def createGraph(pgRNAs, csvout):

    # Open pgRNA target file and extract avg coord distance
    with open(pgRNAs, 'r') as p, open(csvout, 'w') as cA:

        reader = csv.reader(p)
        avg = 0
        count = 0
        cA.write('pgRNA_Location\t' + 'Avg' + '\n')

        for row in reader:

            # if line contains chr#:#-# then handle
            if re.match('chr\S+:\d+-\d+', ' '.join(row)):

                target = ''.join(row)
                
                if avg != 0:

                    avg = avg//count
                    cA.write(target + '\t' + str(avg) + '\n')
                    avg = 0
                    count = 0

            elif re.match('\[', ','.join(row)):

                count += 1

                avg += int(row[6].strip(']'))

    # # calculate pgRNA Distance
    # avg = pd.read_table('avgVal.csv')
    # fif = avg[avg.Avg < 15000]
    # ten = avg[avg.Avg < 10000]
    # ten.hist(bins=1000)
    # plt.title('Avg < 15kb')
    # plt.xlabel('Base Pairs')
    # plt.ylabel('Number of pgRNAs')
    # plt.show()
        
# parse pgRNA.txt and copy pairs used in analysis
def createList(pgRNAs, listout):

    # open pgRNA targets used and extract pairs
    with open(pgRNAs, 'r') as p, open(listout, 'w') as lO, open('pgRNA.txt', 'r') as pg:

        reader = csv.reader(pg)
        pgRNA = csv.reader(p)

        for row in reader:

            if re.match('chr\S+:\d+-\d+', ' '.join(row)):

                for val in pgRNA:
                    
                    if row in val:

                        match = True

                if target in pgRNA:

                    print(target)

## merge lncRNAs used for pgRNA creation and identify non-merged lncRNAs
## three files will be used, 'pgRNA.txt' (contains top 5 pgRNAs per-chr location
## 'jennie_new_wout_bidpl.gtf' (lncRNAs used for Jennie CRISPR deletion)
## and 'crispra_screen.csv' (gRNAs from CRISPRa screen paper)
def mergeId(pgRNAs, gtf, crispra):
    
    # extract chromosome information and save to list
    # this will be used to identify pgRNAs actually found
    with open('pgRNA.txt', 'r') as p, open('lncRNAs_used.csv', 'w') as l:
        for line in p:
            if re.match('chr\S+:\d+-\d+', line):
                l.write(line)

    # extract chromosome and name information and save to file
    df = gtf.dataframe('jennie_new_w_bidpl.gtf')
    df_jen = df['seqname'] + ':' + df['start'] + '-' + df['end'] + '\t' + df['gene_id'] + '\t' + df['gene_name']
    df_jen.to_csv('jennie_w_gtf.csv', index=False)

    # merge 'pgRNA' chromosome location information with chromosome and name from jennie gtf file
    pgRNA = pd.read_table('lncRNAs_used.csv', names=['chrloc'])
    jenGTF = pd.read_table('jennie_w_gtf.csv',  names=['chrloc', 'Ensembl', 'Alias'])
    pg_found = pd.merge(pgRNA, jenGTF, how='inner', on='chrloc')
    pg_found = pg_found.set_index(pg_found[pg_found.columns[2]].str.split('.', expand=True)[0])[['chrloc', 'Ensembl']]

    # extract 'linc_name' from 'crispra_screen.csv' and set to dataframe
    crispra = pd.read_table('crispra_screen.csv', sep=',')
    crispN = crispra[['linc_name']].rename(columns={'linc_name': 'Alias'})
    crispN = crispN.set_index(crispN[crispN.columns[0]].str.split('.', expand=True)[0])
    crispN = crispN[~crispN.index.duplicated(keep='first')]
    merge_new = pd.merge(pg_found, crispN, how='inner', left_index=True, right_index=True)

## Create merged pgRNA file with chr location, id, alias and both gRNAs with names
def finalList(pgRNAs, gtf, final):

    # extract chromosome information and save to list
    # this will be used to identify pgRNAs actually found    
    with open(pgRNA, 'r') as p, open('lncRNAs_used.csv', 'w') as l:
        for line in p:
            if re.match('chr\S+:\d+-\d+', line):
                l.write(line)

    # extract chromosome and name information and save to file
    df = gtf.dataframe('gencode_v22_lncRNA_transcripts.gtf')
    df_jen = df['seqname'] + ':' + df['start'] + '-' + df['end'] + '\t' + df['transcript_id'] + '\t' + df['gene_name']
    df_jen.to_csv('jennie_w_gtf.csv', index=False)

    # merge 'pgRNA' chromosome location information with chromosome and name from gtf file
    pgRNA = pd.read_table('lncRNAs_used.csv', names=['chrloc'])
    jenGTF = pd.read_table('jennie_w_gtf.csv', names=['chrloc', 'Ensembl', 'Alias'])
    pg_found = pd.merge(pgRNA, jenGTF, how='inner', on='chrloc')
    pg_found.to_csv('final_lncRNA_names.tsv', sep='\t', index=False)

    # merge chromosome location information onto pgRNA file
    pg_merge = pd.merge(pgRNA, pg_found, how='outer', on='chrloc')
    pg_merge.to_csv('pgRNA_merged.tsv', sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)

    # parse 'pgRNA_merged.tsv' and reorder information to Jennie's specifications
    with open('pgRNA_merged.tsv', 'r') as pm, open(final, 'w') as f:
        for pm_line in pm:

            # if line contains chromosome coordinate information
            # then safe information to variable
            if re.match('chr\S+:\d+-\d+', pm_line):
                 pm_loc = pm_line
                 mark = False
                 val = 0
            elif mark == False:
                f.write(pm_loc.strip() + '\t' + '\t'.join(pm_line.split(',')[0:4]).replace('["', '').replace('"', '') + '\n')
                val += 1
                if val > 5:
                    mark = True
            elif mark == True:
                pass
                    
    
## script to identify if the previous sequence end site is greater than
## the current sequence start site.  If that is the case then the transcripts
## are said to be overlapping otherwise
def idOverlap(pgRNA, outfile):

    # parse pgRNA file and identify start and end sites
    with open(pgRNA, 'r') as p, open(outfile, 'w') as o:
        prev_line = False
        compare_ready = False
        total = 0
        
        for line in p:

            # split current line if its the chr coordinates and set to prev_line
            if re.match('chr\S+:\d+-\d+', line) and prev_line == False:

                prev_line = line.replace(':', '-').strip().split('-')

            # split current line if its the chr coordinates and prev_line has a value
            elif re.match('chr\S+:\d+-\d+', line) and prev_line != False:

                current_line = line.replace(':', '-').strip().split('-')
                compare_ready = True

            elif compare_ready == True:

                if current_line[0] == prev_line[0] and current_line[1] < prev_line[2]:

                    total +=1

                    #print(str(current_line) + ' < ' + str(prev_line))
                    print(current_line[1] + ' < ' + prev_line[2])
                    o.write(str(current_line) + ' < ' + str(prev_line) + '\n\n')

                compare_ready = False
                prev_line = current_line
        print('Total Problem Sites: ' + str(total))
            
if __name__ == "__main__":

    pgRNA_clean('GuideScan_batch_output.csv', 'top5_test.csv')
    createGraph('pgRNA.txt', 'avgVal.csv')  
        
