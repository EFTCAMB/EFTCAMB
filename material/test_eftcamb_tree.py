###############################################################################
# import modules:
###############################################################################

import sys, platform, os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import integrate
import copy

###############################################################################
# import CAMB:
###############################################################################

#here = os.path.dirname(os.path.abspath(__file__))
here = './'
camb_path = os.path.realpath(os.path.join(os.getcwd(),here))
sys.path.insert(0,camb_path)
import camb
camb.set_feedback_level(10)
from camb import model, initialpower
from camb.baseconfig import CAMBError

#Function to check parameters (it is a prototype still ... but it gets all the model with a double loop)

def check_params(eftcamb_params):

    temp = {}
    flag = ['stability','time','level','priors']

    check = False

    modelflag = ['PureEFTmodel','AltParEFTmodel','DesignerEFTmodel','FullMappingEFTmodel']

    while True :

        # I couldn't find a better way than these flags to come out from the loops (may be improved)

        if temp and temp.get('DesignerEFTmodel',0) == 0 :
            for KEY in temp:
                if ( 'model' in KEY) and temp[KEY] >=4:
                    check = True

        if (temp and (temp.get('DesignerEFTmodel',0) != 0 or temp.get('FullMappingEFTmodel',0) != 0)):
            if temp.get('FullMappingEFTmodel',0) == 4:
                if temp['potential_model'] != 0 :
                    check = True
            else :
               check = True
        if not temp :
          temp =copy.deepcopy(eftcamb_params)

        if check :
            print(temp)
            return temp
            break

        pars = camb.set_params(H0=67.3,**temp)
        read_par = pars.EFTCAMB.read_parameters()

        # This for loop controls the filling of the temp dictionary and assure that only parameter flag are updated
        # in the loop.

        #print(read_par)

        for key in read_par:
            if key not in flag :
                    if (key not in  ['EFTwDE','EFTflag','RPHwDE','potential_model']) and type(read_par[key]) != bool:
                        if key not in modelflag:
                        #print(key,read_par[key],read_par[key] == 0 and 'EFTwDE' != key )
                            temp.update({key : read_par[key]+1})
                        else:
                            temp.update({key : read_par[key]})

#Now we just use python error handling to obtain the whole tree

modelflag = ['PureEFTmodel','AltParEFTmodel','DesignerEFTmodel','FullMappingEFTmodel']
_params=[]

for j,i in enumerate(modelflag):
    for xx in range(1,10):
        if i == modelflag[0] and xx == 2:
            break
        if i == modelflag[-1] and xx == 4:
            eftcamb_params = {'EFTflag':j+1,i:xx,'potential_model': 0}
            for kk in range(1,10):
                eftcamb_params['potential_model'] = kk
                try :
                    params=check_params(eftcamb_params = eftcamb_params)
                    _params.append(list(params.keys()))
                except ValueError:
                    #print(params)
                    break
        else :
            eftcamb_params = {'EFTflag':j+1,i:xx}
            try :
                params=check_params(eftcamb_params = eftcamb_params)
                _params.append(list(params.keys()))
            except ValueError:
                #print(params)
                break
