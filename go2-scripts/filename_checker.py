#!/usr/bin/env python3

# Dependencies
try:
    import os
    import sys
    import re
    import argparse
    import pandas as pd
    import datetime as dt
    from pathlib import Path
except ImportError:
    print('Some modules not installed')
    print("Check and install modules: os, sys, re, argparse, pandas,\
            datetime, pathlib")
    sys.exit('ImportError')

#----------------------------------------------------------------------------------------
#     ***************           Functions                ****************
#========================================================================================
def get_fields(inFile):
    """
    Splits file name from path and parses filename into words.
    returns a dict or list
    input: filename and pathe
           text
    output: boolen & dict
            turple
    """
    # Expected filename naming convention
    convention1= ['study', 'pheno', 'ethnicity', 'ncases', 'ncontrols','software',\
                  'panel','build','date','analyst']
    convention2= ['study', 'pheno', 'ethnicity', 'ncases', 'ncontrols','software',\
                  'panel','build','date','analyst','gc','ff']
    #=================================================================================

    # Split path and filename
    (npath, filename) = os.path.split(inFile) 

    # Split filename
    splitname= filename.split('_')
    splitResult= splitname[:-1]
    splitResult.append(splitname[-1].split('.')[0])
    nfields= len(splitResult)

    #Get absolute path of file
    path=str(Path(inFile).absolute())

    outDat={'path':path}
    if nfields==10:
        nameConvention=True # filename valid
        for f, name in zip(convention1, splitResult):
            outDat[f]= name
        outDat['gc']= 'EMPTY'
        outDat['ff']= 'EMPTY'
        return nameConvention, outDat

    elif nfields==12:
        nameConvention=True
        for f, name in zip(convention2, splitResult):
            outDat[f]= name

        return nameConvention, outDat
    else:
        nameConvention=False
        print('\n\n\t******* Parsing filename *******')
        print(f"ERROR: Number of fields {{nfields}} in file name must be 10 or 12")
        print(f"Number of fields _nfields_ in file name is {nfields}")
        print(f'File name :{filename}')

        sys.exit(ValueError)

    return nameConvention, splitResult

# Fxn Call
#------------------------------
#t1=get_fields(inFile=path1)
###########################################################################################
def check_pheno(infield, phenlist, prefixList):
    prefix= [x.split(',') for x in prefixList][0]
    prefix= [x.replace(' ', '') for x in prefix]

    ##Get all combinations of prefix

    allpheno= [x.split(',') for x in phenlist][0]
    allpheno= [x.replace(' ', '') for x in allpheno]

    outField= 'None'
    for p in prefix:
        if (infield.startswith(p)):
            new_phenlist= [p+x for x in allpheno]
            partn= rf"\b({infield})\b"
            new_phenoAll= ' '.join(new_phenlist)
            check=re.search(partn, new_phenoAll, flags=re.IGNORECASE)
            if check:
                outField= infield
                break

    if outField != infield:
        outField='ERROR'

        print(f'Warning: {key_name} {{infield}} _{infield}_ not found')

        print(f'Warning: pheno  {{infield}} _{infield}_ not found')

    return outField
# FUN Call
#test=check_pheno(infield, phenoList, prefixList)
#-----------------------------------------------------------------------------------

def check_fields(nameFields, conFile=None):

    """
    Verifies that file name meets naming convention
    Returns a dict when filename is ok, and warning otherwise

    nameFields:  dict
                 from get_fields fxn
    output: dict
    """

    # Get user values
    def read_vals(conFile):
        with open(conFile) as f:
            lines = [line.rstrip() for line in f]
        return lines


    # List of field elements expected to be present in file names
    # Can be updated or alternatively set as argument
    #Set Authorised Values

    if conFile is not None:
        # Get user values
        lines=read_vals(conFile)

        # Gett List values
        for line in lines:
            list_vals= line.split(':')[1:]

            if line.startswith('pheno'):
                phenoList=list_vals

            elif line.startswith('ethnicity'):
                ethnicList= list_vals

            elif line.startswith('software'):
                softwareList= list_vals

            elif line.startswith('panel'):
                panelList= list_vals

            elif line.startswith('build'):
                buildList= list_vals
            elif line.startswith('gc'):
                gcList= list_vals
            elif line.startswith('ff'):
                ffList=list_vals
            elif line.startswith('prefix'):
                prefixList=list_vals
    else:
        phenoList= ['ALLOA', 'HIP', 'MALEKNEE', 'HIPKNEE', 'HAND', 'SPINE', 'THUMB',\
            'MALESPINE', 'FINGER', 'TJR', 'THR', 'TKR']
        ethnicList= ['EA', 'AA', 'HL', 'EAS', 'SAS', 'AU']
        softwareList=['SNPTEST','BOLTLMM','SAIGE','GEMMA','PLINK','RVTEST','OTHER']
        panelList=['BESPOKE','HRC','1000G','NA']
        buildList=[36, 37, 38 ]
        gcList=['GCYes', 'GCNo']
        ffList=['FFYes', 'FFNo']
        prefixList=['MALE', 'FEMALE', 'EARLY', 'EARLYMALE', 'EARLYFEMALE']
    #========================================================================

    outCheck={}
    for field in nameFields:
        key_name= field.upper()
        infield= nameFields[field]
        #print(field.upper())

        key_name= field.upper()
        infield= nameFields[field]
        
        
        if field == 'pheno':
            partn= rf"\b({infield})\b"
            phenoAll= ' '.join(phenoList)
            check=re.search(partn, phenoAll, flags=re.IGNORECASE)
            if check:
                outField= infield
            else:
                ## Check pheno with combinations
                # Run check_pheno()
                outField=check_pheno(infield, phenoList, prefixList)

        elif field =='ethnicity':
            partn= rf"\b({infield})\b"
            ethnicityAll= ' '.join(ethnicList)
            check=re.search(partn, ethnicityAll, flags=re.IGNORECASE)
            if check:
                outField=infield
            else:
                outField='ERROR'
                print(f'Warning: {key_name} {{infield}} _{infield}_ not found')

        elif field=='ncases':
            check=re.search(r"(CASE)S?$", infield, flags=re.IGNORECASE)

            if check:
                # Where field ends with suffix (Cases)
                check_num= re.findall(r"\d+", infield)[0]

                try:
                    cases= int(check_num)
                    outField=  str(cases)
                except ValueError:
                    outField= 'ERROR'
                    print(f'Warning: {key_name} {{infield}} _{infield}_ is invalid')
            else:
                # Where field has no suffix (Cases)
                try:
                    cases= int(infield)
                    outField=  infield
                except ValueError:
                    outField= 'ERROR'
                    print(f'Warning: {key_name} {{infield}} _{infield}_ is invalid')

        elif field=='ncontrols':
            check=re.search(r"(CONTROL)S?$", infield, flags=re.IGNORECASE)

            if check:
                # Where field ends with suffix (Controls)
                check_num= re.findall(r"\d+", infield)[0]

                try:
                    controls= int(check_num)
                    outField=  str(controls)
                except ValueError:
                    outField= 'ERROR'
                    print(f'Warning: {key_name} {{infield}} _{infield}_ is invalid')
            else:
                # Where field has no suffix (Controls)
                try:
                    controls= int(infield)
                    outField=  infield
                except ValueError:
                    outField= 'ERROR'
                    print(f'Warning: {key_name} {{infield}} _{infield}_ is invalid')

        elif field=='software':
            partn= rf"\b({infield})\b"
            softwareAll= ' '.join(softwareList)
            check= re.search(partn, softwareAll, flags=re.IGNORECASE)
            if check:
                outField= infield
            else:
                outField='ERROR'
                print(f'Warning: {key_name} {{infield}} _{infield}_ not found')

        elif field== 'panel':
            partn= rf"\b({infield})\b"
            panelAll= ' '.join(panelList)
            check= re.search(partn, panelAll, flags=re.IGNORECASE)
            if check:
                outField= infield
            else:
                outField='ERROR'
                print(f"Warning: {key_name} {{infield}} _{infield}_ not found")

        elif field=='build':
            partn= rf"\b({infield})\b"
            temp=[str(x) for x in buildList]
            buildAll= ' '.join(temp)
            check= re.search(partn, buildAll, flags=re.IGNORECASE)
            if check:
                outField=infield
            else:
                outField='ERROR'
                print(f'Warning: {key_name} {{infield}} _{infield}_ not found')

        elif field== 'date':
            try:
                if type(dt.datetime.strptime(infield, "%d%m%y").date()) is dt.date:
                    outField= infield
                    print('Date True1: DDMMYY')

            except ValueError:
                #print('False1')
                pass # check month-date-year format
                try:
                    if type(dt.datetime.strptime(infield, "%m%d%y").date()) is dt.date:
                        outField= infield
                        print('Date True2: MMDDYY')
                except ValueError:
                    #print('False2')
                    pass # check DDMMYY format
                    try:
                        if type(dt.datetime.strptime(infield, "%d%m%Y").date()) is dt.date:
                            outField= infield
                            print('Date True3: DDMMYYYY')
                    except ValueError:
                        #print('False3')
                        pass
                        try:
                            if type(dt.datetime.strptime(infield, "%m%d%Y").date()) is dt.date:
                                outField= infield
                                print('Date True4: MMDDYYYY')
                        except ValueError:
                            #print('False4')
                            outField= 'ERROR'
                            print(f'Warning: {key_name} {{infield}} _{infield}_ format is invalid')
                            print('Date field invalid. Date format is neither DDMMYY or MMDDYY or DDMMYYYY or MMDDYYYY')

        elif field=='analyst':
            outField= infield
        elif (field=='gc') & (infield !='EMPTY'):
            partn= rf"\b({infield})\b"
            gcAll= ' '.join(gcList)
            check= re.findall(partn, gcAll, flags=re.IGNORECASE)
            if len(check) > 0:
                outField= check[0]
            else:
                outField= 'ERROR'
                print(f'Warning: {key_name} {{infield}} _{infield}_ not found')

        elif (field=='gc') & (infield =='EMPTY'):
            # Follows first naming convention with 10 fields
            outField= infield

        elif (field=='ff') & (infield !='EMPTY'):
            partn= rf"\b({infield})\b"
            ffAll= ' '.join(ffList)
            check= re.findall(partn, ffAll, flags=re.IGNORECASE)
            if len(check) > 0:
                outField= check[0]
            else:
                outField= 'ERROR'
                print(f'Warning: {key_name}(Fudge Factor) {{infield}} _{infield}_ not found')

        elif (field=='ff') & (infield =='EMPTY'):
            # Follows naming convention with 10 fields
            outField= infield

        elif field =='study':
            outField= infield
        elif field== 'path':
            outField= infield


        outCheck[key_name]= outField
        #print(outField)


    if re.search(r'\b(ERROR)\b', ' '.join(outCheck.values()), flags=re.I):
        outCheck['ISERROR'] ='YES'
    else:
        outCheck['ISERROR'] = 'NO'

    return outCheck

#-----------------------------------------------------------------------------------------------------
#######################################################################################################
#== input args  ====
#inFile= 'args1'
#conFile= 'args2'
#outputFile= 'args3'

#usage
# ./filename_checker [FILEPATHS] [-v FILEVALS] [-o OUTNAME]
#-------------------------------------------------------------------------------------------------------

### Parse user input arguments
parser= argparse.ArgumentParser()
parser.add_argument("FILEPATHS", help="text file with list of file paths.")
parser.add_argument("-v", "--filevals", help="text file with authorised values for file name fields.\
                                            File name submitted as string")
parser.add_argument("-o", "--outname", help="string text to be used as output file name")
args= parser.parse_args()


# args 1
if args.filevals:
        confile= True
        print('Config File:', args.filevals)
else: confile= False

# args 2
try:
    if re.search(r'\b(.list.txt)\b', args.FILEPATHS, flags=re.I):
        inlist= True
    else:
        inlist= False
except TypeError:
    print('ERROR: List of file paths not found')
    sys.exit(TypeError)

# args 3
if args.outname:
    outputFile= args.outname
else:
    outputFile='files.checkout'


#----------------------------------------------------------------------------------------------------------

### Main
# Read FILEPATHSs into a list
with open(args.FILEPATHS) as f:
    lines = [line.rstrip() for line in f]

# Check each FILEPATHS
outAll= pd.DataFrame()
for line in lines:
    ### 1. Get fields
    outfields=get_fields(line)

    (head, tail) = os.path.split(line)
    print('\n\n\t******* Parsing filename *******')

    ### 2. Check fields
    if (outfields[0]==1) & (confile==1):
        outcheck= check_fields(nameFields=outfields[1], conFile=args.filevals)
        outframe=pd.DataFrame.from_dict([outcheck])
        if re.search(r'\b(ERROR)\b', ' '.join(outcheck.values()), flags=re.I):
            print('\n\tNote: File name not valid')
        else:   print('\tFile name is valid')
        print(f'File name: {tail}')

    elif outfields[0]==1 : # No. of fields in name as expected and no config file
        outcheck= check_fields(nameFields=outfields[1])
        outframe=pd.DataFrame.from_dict([outcheck])
        if re.search(r'\b(ERROR)\b', ' '.join(outcheck.values()), flags=re.I):
            print('\tNote: File name not valid')
        else: print('\tFile name is valid')
        print(f'File name: {tail}')

    outAll= outAll.append(outframe)

### 3. Save output
outAll.to_csv(f"{outputFile}.csv", index=False)
print(f"Output file: {outputFile}.csv")

