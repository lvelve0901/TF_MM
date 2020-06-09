# Program to check if we can see if the base pair is in a helix
# First, take a control that is within a helix
import os
import pandas as pd
import time
mode = 0

def checkinteger(test_char):
    ''' Check if test_char represents an integer '''
    try:
        int(test_char)
        return True
    except ValueError:
        return False

def filter_helix(input_file,output_file):
    global mode
    # Mode - Search for canonical duplexes
    mode = 1
    # Search for mispaired structures
    #mode = -1

    if mode == 1:
        print "Searching for canonical dupex contexts "
        time.sleep(0.5)
    elif mode == -1:
        print "Searching for mispaired duplex contexts "
        time.sleep(0.5)
    
    base_pairs = ["AT", "aT", "At", "at", "TA", "tA", "Ta", "ta", "GC", "gC", "Gc", "gc", "CG", "cG", "Cg", "cg", "AU", "Au", "aU", "au", "UA", "uA", "Ua", "ua"]
    input_bp = pd.read_csv('%s'%input_file)
    all_bp = pd.read_csv('../../PairTable_nmr.csv')
    status_bp = []
 
    # Convert the chain data types to str (for both chains)
    all_bp['chid_1'] = all_bp['chid_1'].astype('str')
    all_bp['chid_2'] = all_bp['chid_2'].astype('str')
    all_bp['resi_1'] = all_bp['resi_1'].astype('str')
    all_bp['resi_2'] = all_bp['resi_2'].astype('str')
    input_bp['chid_1'] = input_bp['chid_1'].astype('str')
    input_bp['chid_2'] = input_bp['chid_2'].astype('str')
    input_bp['resi_1'] = input_bp['resi_1'].astype('str')
    input_bp['resi_2'] = input_bp['resi_2'].astype('str')

    #for dummy in range(len(all_bp.chid_1)):
    #    all_bp.chid_1[dummy] = str(all_bp.chid_1[dummy])
    #for dummy in range(len(all_bp.chid_2)):
    #    all_bp.chid_1[dummy] = str(all_bp.chid_2[dummy])
    print "Done"
    #all_bp.loc[:, 'chid_1'] = str(all_bp.loc[:, 'chid_1'])
    #all_bp.loc[:, 'chid_2'] = str(all_bp.loc[:, 'chid_2'])

    # For each CC base pair
    for ag_basepair in range(len(input_bp)):
    #for ag_basepair in range(157, 160):
        # Search in the list of all base pairs for the corresponding pdb  
        # Fetch all the base pairs in that pdb
        pdb_bp = all_bp.loc[all_bp.pdbid.isin([input_bp.pdbid[ag_basepair]])]
        #print "CHECK ", len(pdb_bp)
        #print pdb_bp
        # Search over the above list for instances in which 
        #print input_bp.chid_1[ag_basepair], input_bp.resi_1[ag_basepair], input_bp.chid_2[ag_basepair], input_bp.resi_2[ag_basepair]
        #print type(input_bp.chid_1[ag_basepair]), type(input_bp.resi_1[ag_basepair]), type(input_bp.chid_2[ag_basepair]), type(input_bp.resi_2[ag_basepair])
        # Check for data types
        if ((checkinteger(input_bp.resi_1[ag_basepair]) == False) or (checkinteger(input_bp.resi_2[ag_basepair]) == False)):
            print input_bp.pdbid[ag_basepair], "Check Manually"
            status_bp.append("Check Manually")
            continue
        #print input_bp.pdbid[ag_basepair], type(input_bp.chid_1[ag_basepair]), type(input_bp.chid_2[ag_basepair])
    
        # Now check for one base pair above 
        bp_above = pdb_bp.loc[(((pdb_bp.chid_1 == input_bp.chid_1[ag_basepair]) & (pdb_bp.chid_2 == input_bp.chid_2[ag_basepair]) & ((pdb_bp.resi_1 == str(int(input_bp.resi_1[ag_basepair]) - 1))|(pdb_bp.resi_1 == int(input_bp.resi_1[ag_basepair]) - 1)) & ((pdb_bp.resi_2 == str(int(input_bp.resi_2[ag_basepair]) + 1)) | (pdb_bp.resi_2 == int(input_bp.resi_2[ag_basepair]) + 1))) | ((pdb_bp.chid_1 == input_bp.chid_2[ag_basepair]) & (pdb_bp.chid_2 == input_bp.chid_1[ag_basepair]) & ((pdb_bp.resi_1 == str(int(input_bp.resi_2[ag_basepair]) + 1)) | (pdb_bp.resi_1 == int(input_bp.resi_2[ag_basepair]) + 1)) & ((pdb_bp.resi_2 == str(int(input_bp.resi_1[ag_basepair]) - 1)) | (pdb_bp.resi_2 == int(input_bp.resi_1[ag_basepair]) - 1))))]
        #print input_bp.chid_1[ag_basepair], input_bp.chid_2[ag_basepair], int(input_bp.resi_1[ag_basepair]) - 1, int(input_bp.resi_2[ag_basepair]) + 1, input_bp.chid_2[ag_basepair], input_bp.chid_1[ag_basepair], int(input_bp.resi_2[ag_basepair]) + 1, int(input_bp.resi_1[ag_basepair]) - 1
    
        # Now check for one base pair beloe
        bp_below = pdb_bp.loc[((((pdb_bp.chid_1 == input_bp.chid_1[ag_basepair]) & (pdb_bp.chid_2 == input_bp.chid_2[ag_basepair]) & ((pdb_bp.resi_1 == str(int(input_bp.resi_1[ag_basepair]) + 1)) | (pdb_bp.resi_1 == int(input_bp.resi_1[ag_basepair]) + 1)) & ((pdb_bp.resi_2 == str(int(input_bp.resi_2[ag_basepair]) - 1)) | (pdb_bp.resi_2 == int(input_bp.resi_2[ag_basepair]) - 1)))) | ((pdb_bp.chid_1 == input_bp.chid_2[ag_basepair]) & (pdb_bp.chid_2 == input_bp.chid_1[ag_basepair]) & ((pdb_bp.resi_1 == str(int(input_bp.resi_2[ag_basepair]) - 1)) | (pdb_bp.resi_1 == int(input_bp.resi_2[ag_basepair]) - 1)) & ((pdb_bp.resi_2 == str(int(input_bp.resi_1[ag_basepair]) + 1)) | (pdb_bp.resi_2 == int(input_bp.resi_1[ag_basepair]) + 1))))]

    
        # Now check for two base pairs above
        bp_two_above = pdb_bp.loc[(((pdb_bp.chid_1 == input_bp.chid_1[ag_basepair]) & (pdb_bp.chid_2 == input_bp.chid_2[ag_basepair]) & ((pdb_bp.resi_1 == str(int(input_bp.resi_1[ag_basepair]) - 2)) | (pdb_bp.resi_1 == int(input_bp.resi_1[ag_basepair]) - 2)) & ((pdb_bp.resi_2 == str(int(input_bp.resi_2[ag_basepair]) + 2)) | (pdb_bp.resi_2 == int(input_bp.resi_2[ag_basepair]) + 2))) | ((pdb_bp.chid_1 == input_bp.chid_2[ag_basepair]) & (pdb_bp.chid_2 == input_bp.chid_1[ag_basepair]) & ((pdb_bp.resi_1 == str(int(input_bp.resi_2[ag_basepair]) + 2)) | (pdb_bp.resi_1 == int(input_bp.resi_2[ag_basepair]) + 2)) & ((pdb_bp.resi_2 == str(int(input_bp.resi_1[ag_basepair]) - 2)) | (pdb_bp.resi_2 == int(input_bp.resi_1[ag_basepair]) - 2))))]
        
    
        # Now check for two base pairs below
        bp_two_below = pdb_bp.loc[(((pdb_bp.chid_1 == input_bp.chid_1[ag_basepair]) & (pdb_bp.chid_2 == input_bp.chid_2[ag_basepair]) & ((pdb_bp.resi_1 == str(int(input_bp.resi_1[ag_basepair]) + 2)) | (pdb_bp.resi_1 == int(input_bp.resi_1[ag_basepair]) + 2)) & ((pdb_bp.resi_2 == str(int(input_bp.resi_2[ag_basepair]) - 2)) | (pdb_bp.resi_2 == int(input_bp.resi_2[ag_basepair]) - 2))) | ((pdb_bp.chid_1 == input_bp.chid_2[ag_basepair]) & (pdb_bp.chid_2 == input_bp.chid_1[ag_basepair]) & ((pdb_bp.resi_1 == str(int(input_bp.resi_2[ag_basepair]) - 2)) | (pdb_bp.resi_1 == int(input_bp.resi_2[ag_basepair]) - 2)) & ((pdb_bp.resi_2 == str(int(input_bp.resi_1[ag_basepair]) + 2)) | (pdb_bp.resi_2 == int(input_bp.resi_1[ag_basepair]) + 2))))]

     
        #print
        #print input_bp.pdbid[ag_basepair]
        #print bp_above, len(bp_above)
        #print bp_below, len(bp_below)
        #print bp_two_above, len(bp_two_above)
        #print bp_two_below, len(bp_two_below)
    
        # The > sign is for handling symmetric constructs
        if (mode == 1):
            if (((len(bp_above) >= 1) and (len(bp_below) >= 1) and (len(bp_two_above) >= 1) and (len(bp_two_below) >= 1))):
                # Get all the base pair names
                names_nbors = []
                for name in bp_above.bp_name:
                    names_nbors.append(name)
                for name in bp_below.bp_name:
                    names_nbors.append(name)
                for name in bp_two_above.bp_name:
                    names_nbors.append(name)
                for name in bp_two_below.bp_name:
                    names_nbors.append(name)

                # Check if they are non-canonical
                mismatch_status = 0
                for name in names_nbors:
                    if name not in base_pairs:
                        mismatch_status = 1

                # Write o/p accordingly
                if mismatch_status == 1:
                    print input_bp.pdbid[ag_basepair], "Mismatched Duplex"
                    status_bp.append("Mismatched Duplex")
                else:
                    print input_bp.pdbid[ag_basepair], "Duplex"
                    status_bp.append("Duplex")
            else:
                print input_bp.pdbid[ag_basepair], "Not categorizable"
                status_bp.append("Not categorizable")

        elif mode == -1:
            if (((len(bp_above) >= 1) and (len(bp_below) >= 1) and (len(bp_two_above) >= 1) and (len(bp_two_below) >= 1))):
                # Get all the base pair names
                names_nbors = []
                names_nbors_far_above = []
                names_nbors_far_below = []
                names_nbors_close = []
                for name in bp_above.bp_name:
                    names_nbors.append(name)
                    names_nbors_close.append(name)
                for name in bp_below.bp_name:
                    names_nbors.append(name)
                    names_nbors_close.append(name)
                for name in bp_two_above.bp_name:
                    names_nbors.append(name)
                    names_nbors_far_above.append(name)
                for name in bp_two_below.bp_name:
                    names_nbors.append(name)
                    names_nbors_far_below.append(name)

                canonical_status = 0
                canonical_status_close = 0
                canonical_status_far_above = 0
                canonical_status_far_below = 0
                for name in names_nbors:
                    if name in base_pairs:
                        canonical_status = 1

                for name in names_nbors_close:
                    if name in base_pairs:
                        canonical_status_close = 1

                for name in names_nbors_far_above:
                    if name in base_pairs:
                        canonical_status_far_above = 1

                for name in names_nbors_far_below:
                    if name in base_pairs:
                        canonical_status_far_below = 1

                # Write o/p accordingly
                if (canonical_status == 0) and (canonical_status_close == 0) and (canonical_status_far_above == 0) and (canonical_status_far_below == 0):
                    print input_bp.pdbid[ag_basepair], "5 Mispairs"
                    status_bp.append("5 Mispairs")
                elif (canonical_status == 1) and (canonical_status_close == 0) and ((canonical_status_far_above == 0) or (canonical_status_far_below == 0)):
                    print input_bp.pdbid[ag_basepair], "4 Mispairs"
                    status_bp.append("4 Mispairs")
                elif (canonical_status == 1) and (canonical_status_close == 0):
                    print input_bp.pdbid[ag_basepair], "3 Mispairs"
                    status_bp.append("3 Mispairs")
                else:
                    print input_bp.pdbid[ag_basepair], "Canonical Duplex"
                    status_bp.append("Canonical Duplex")
            else:
                print input_bp.pdbid[ag_basepair], "Not categorizable"
                status_bp.append("Not categorizable")
 
                   

    print
    # Once done
    print len(status_bp)
    print len(input_bp.pdbid)
    input_bp.insert(1, 'status_bp', status_bp)
    input_bp = input_bp.loc[input_bp['status_bp'] == 'Duplex']
    input_bp.to_csv('%s'%output_file,index=False)

def main():
    
    filelist = ["DNA_mismatch_GG_nmr"]
    #filelist = ["DNA_mismatch_GT_nmr","DNA_mismatch_AC_nmr","DNA_mismatch_GG_nmr","DNA_mismatch_AG_nmr","DNA_mismatch_AA_nmr","DNA_mismatch_TT_nmr","DNA_mismatch_CT_nmr","DNA_mismatch_CC_nmr"]

    for filename in filelist:
        input_file = "%s.csv"%filename
        output_file = "%s_automatic.csv"%filename
        filter_helix(input_file,output_file)
        print "Finished: %s" %filename


if __name__ == "__main__":
    main()

