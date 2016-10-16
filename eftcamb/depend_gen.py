#----------------------------------------------------------------------------------------
#
# This file is part of EFTCAMB.
#
# Copyright (C) 2013-2016 by the EFTCAMB authors
#
# The EFTCAMB code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file eftcamb/LICENSE at
# the top level of the EFTCAMB distribution.
#
#----------------------------------------------------------------------------------------

#
# This file contains a python script that automatically generates the source file
# dependencies for the EFTCAMB code.
# This heavily relies on a python script by David Dickinson and Peter Hill.
#

#!/usr/bin/python
import os
import re

# hard coded options:
output     = 'eftcamb.dep'
path_out   = '$(OUTPUT_DIR)/'
path_dll   = '$(DLL_DIR)/'

#Definitions
def run( verbose=True, overwrite=None, outdir=None ):

    # some feedback:
    print "\033[032m"+'Generating EFTCAMB dependencies'+"\033[039m"
    
    # grab the EFTCAMB/CAMB paths:
    eftcamb_dir = os.path.dirname(os.path.abspath(__file__))
    camb_dir    = os.path.abspath(os.path.join( eftcamb_dir, os.pardir ))

    # get the EFTCAMB/CAMB files:
    eftcamb_files = []
    for root, dirs, files in os.walk(eftcamb_dir):
        for file in files:
            if root != eftcamb_dir:
                eftcamb_files.append( os.path.join(os.path.basename(root),file) )
            else:
                eftcamb_files.append( file )
    
    camb_files    = os.listdir(camb_dir)

    # filter them to get only the fortran files:
    eftcamb_files = select_sources(eftcamb_files)
    camb_files    = select_sources(camb_files)

    # append path:
    eftcamb_files_temp = []
    for f in eftcamb_files:
        eftcamb_files_temp.append( os.path.join( eftcamb_dir, f ) )
    eftcamb_files = eftcamb_files_temp

    camb_files_temp = []
    for f in camb_files:
        camb_files_temp.append( os.path.join( camb_dir, f ) )
    camb_files = camb_files_temp

    l_camb    = create_file_objs(camb_files)
    l_eftcamb = create_file_objs(eftcamb_files)

    mod2fil_camb      = file_objs_to_mod_dict( file_objs=l_camb )
    mod2fil_eftcamb   = file_objs_to_mod_dict( file_objs=l_eftcamb )

    # we have to write this because we do not want to mess with CAMB dependencies
    dependencies={}
    # start with EFTCAMB:
    for i in l_eftcamb:
        tmp=[]
        for j in i.uses:
            try:
                tmp.append( mod2fil_eftcamb[j.lower()] )
            except:
                try:
                    tmp.append( mod2fil_camb[j.lower()] )
                except:
                    print "\033[031mError\033[039m module \033[032m"+j+"\033[039m not defined in any files. Skipping..."
        if len(tmp)>0:
            dependencies[i.file_name]=tmp

    for i in l_camb:
        tmp=[]
        for j in i.uses:
            try:
                tmp.append( mod2fil_eftcamb[j.lower()] )
            except:
                pass
        if len(tmp)>0:
            dependencies[i.file_name]=tmp

    if verbose:
        for i in dependencies.keys():
            print "\033[032m"+i+"\033[039m depends on :\033[034m"
            for j in dependencies[i]: print "\t"+j
            print "\033[039m"

    if outdir is not None:
        outname = os.path.join( outdir, output )
    else:
        outname = output

    print "\033[032m"+'Saving dependencies to: '+"\033[039m"+outname

    tmp=write_depend(outfile=outname,dep=dependencies,overwrite=True , path=path_out)
    tmp=write_depend(outfile=outname,dep=dependencies,overwrite=False, path=path_dll)

    return dependencies

def write_depend(outfile="makefile.dep",dep=[],overwrite=False, path=""):
    "Write the dependencies to outfile"
    #Test file doesn't exist
    if not(overwrite):
        opt="n"
    else:
        opt="y"

    #Open file
    if opt == "y":
        f = open(outfile,'w')
    elif opt == "n":
        f = open(outfile,'a')

    for i in dep.keys():
        tmp,fil=os.path.split(i)
        stri="\n"+path+fil.split(".")[0]+".o"+" : "
        for j in dep[i]:
            tmp,fil=os.path.split(j)
            stri=stri+" \\\n\t"+path+fil.split(".")[0]+".o"
        stri=stri+"\n"
        f.write(stri)
    f.close()
    return

def select_sources( files, ext=[".f",".f90",".F90"] ):
    "Return all files ending with any of ext in the input list"
    fil=[]
    for i in ext:
        fil.extend(filter(lambda x: x.endswith(i),files))
    return fil

def get_source(ext=[".f",".f90",".F90"]):
    "Return all files ending with any of ext"
    tmp=os.listdir(".")
    fil=[]
    for i in ext:
        fil.extend(filter(lambda x: x.endswith(i),tmp))
    return fil

def create_file_objs(files=None, macros={}):
    l=[]

    if files is None:
        files = get_source()

    for i in files:

        source_file = file_obj()

        source_file.file_name = os.path.basename( i )
        source_file.path      = os.path.dirname( i )
        source_file.uses      = get_uses(i,macros)
        source_file.contains  = get_contains(i)

        l.append(source_file)

    return l

def get_uses(infile=None, macros={}):
    "Return which modules are used in infile after expanding macros"
    p=re.compile("^\s*use\s*(?P<moduse>\w*)\s*(,)?\s*(only)?\s*(:)?.*?$",re.IGNORECASE).match

    uses=[]

    with open(infile,'r') as f:
        t=f.readlines()

    for i in t:
        tmp=p(i)
        if tmp:
            uses.append(tmp.group('moduse').strip())

    # Remove duplicates
    uniq_mods = list(set(uses))

    for i, mod in enumerate(uniq_mods):
        for k, v in macros.items():
            if re.match(k, mod, re.IGNORECASE):
                uniq_mods[i] = mod.replace(k,v)

    return uniq_mods

def get_contains(infile=None):
    "Return all the modules that are in infile"
    p=re.compile("^\s*module\s*(?P<modname>\w*?)\s*$",re.IGNORECASE).match

    contains=[]

    with open(infile,'r') as f:
        t=f.readlines()

    for i in t:
        tmp=p(i)
        if tmp:
            contains.append(tmp.group('modname').strip())

    # Remove duplicates before returning
    return list(set(contains))

def file_objs_to_mod_dict(file_objs=[]):
    "Turn a list of file_objs in a dictionary, containing which modules depend on which files"
    dic={}
    for i in file_objs:
        for j in i.contains:
            dic[j.lower()]=i.file_name
    return dic

def get_depends(fob=[],m2f=[]):
    deps={}
    for i in fob:
        tmp=[]
        for j in i.uses:
            try:
                tmp.append(m2f[j.lower()])
            except:
                print "\033[031mError\033[039m module \033[032m"+j+"\033[039m not defined in any files. Skipping..."

        deps[i.file_name]=tmp

    return deps

class file_obj:
    def __init__(self):
        self.path=None
        self.file_name=None
        self.uses=None
        self.contains=None
        self.depends_on=None


#Script
if __name__ == "__main__":

    import argparse

    # Add command line arguments
    parser = argparse.ArgumentParser(description='Generate Fortran dependencies for EFTCAMB')
    parser.add_argument('-v','--verbose',action='store_true',help='explain what is done')
    parser.add_argument('-o','--outroot', dest='outroot', type=str,help='output directory for the file')

    # Parse the command line arguments
    args = parser.parse_args()

    run( verbose=args.verbose, outdir=args.outroot )
