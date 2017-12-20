#############################################
# IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import os
from scipy import optimize
try:
	import argparse
	argparser=1
except:
	print "Python module >> argparse << is not installed. Argument parsing is not possible. Script will not work."
	argparser=0
import sys
#############################################
# DEFINITION OF FUNCTIONS

def fcn2min(params, x, data):
    """ observed protein shifts, subtract data"""
    global P0
    global n
    Kd = params['Kd']
    Dmax = params['Dmax']
    model = Dmax/2.0 * ( (x*n+1.+n*Kd/P0) - np.sqrt( (x*n+1.+n*Kd/P0)**2.0 - 4.0*x*n ) )
    return model - data

# Function to calculate the simulated curve (for plotting)
def simcurve(x,Kd,Dmax):
    global P0
    global n
    simcurve =  Dmax/2.0 * ( (x*n+1.+n*Kd/P0) - np.sqrt( (x*n+1.+n*Kd/P0)**2.0 - 4.0*x*n ) )
    print len(simcurve)
    return simcurve

# Function used to create a directory
def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete        - regular file in the way, raise an exception        - parent
          directory(ies) does not exist, make them as well    """ 
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        #print "_mkdir %s" % repr(newdir)
        if tail:
            os.mkdir(newdir)
###########################################################

################################################
# M A I N
################################################

# Read data
print "Type -h for help"
# SPECIFY INPUT PARAMETERS. FIRST FOR THE CASE THAT ARGPARSE IS POSSIBLE (MODULE INSTALLED)
if argparser==1:
	#### Handle the input and output arguments
	parser = argparse.ArgumentParser(description='This script fits Kd values from titration data, assuming a global Kd for multiple residues.')
	parser.add_argument('-d','--data', help='Input file name. Format: First line, first col: "0"; following columns residue numbers. All other lines: first column: ligand concentration relative to protein concentration, e.g. "1" for equal molar amounts. Other columns: chemical-shift differences.)',required=True,default='data.txt')
	parser.add_argument('-conc','--proteinconcentration',help='Protein concentration in molar', required=True,default=100e-6)
	parser.add_argument('-n','--stoichiometry',help='Stoichiometry parameter n', required=True,default=1)
	parser.add_argument('-o','--outputdir',help='Name of output directory', required=True,default='output')
	parser.add_argument('-Kd','--Kd_estimate',help='estimated Kd', required=False,default=10e-6)
	parser.add_argument('-Dmax','--Dmaxfile',help='Optional file with Dmax values to be fixed', required=False)
	parser.add_argument('-fixKd','--fixedKdvalue',help='A uniform Kd value to be fixed when fitting each residue', required=False)
	args = parser.parse_args()
	outputdir=args.outputdir
	shiftfile=args.data
	P0=float(args.proteinconcentration)
	n=int(args.stoichiometry)
	Kd_init=float(args.Kd_estimate)
	Dmaxfile=args.Dmaxfile
	fixedKd=args.fixedKdvalue

_mkdir(outputdir)
#########################################
# Read the input file containing the chemical-shift data	
data=np.loadtxt(shiftfile,skiprows=1)
data = data.transpose()
x=data[0]
data=data[1:]
# Read the first line of the file, which contains the residues numbers
residues=np.loadtxt(shiftfile)[0]
residues=residues[1:]

########################################
# Optionally read a file that contains the Dmax values
try:
	Dmaxlist=np.loadtxt(Dmaxfile)
	print "Dmax values will be fixed to values from file ",Dmaxfile
	fixDmax=1	# flag to fix the Dmax values
except:
	print "No file with Dmax values specified. Dmax values will be fitted."
	fixDmax=0

#########################################
# Optionally fix the Kd
if fixedKd:
	fixedKd=float(fixedKd)

#########################################
# Auxiliary stuff for writing/plotting data
simpoints=100
x_sim=np.linspace(min(x),max(x)*2,simpoints)
simdata=np.zeros((len(residues)+1,simpoints+1))
simdata[0][1:]=x_sim
simdata=np.transpose(simdata)
simdata[0][1:]=residues
simdata=np.transpose(simdata)

###########################################################################
# FITS
###########################################################################
summaryfile = open(str(outputdir)+'/summaryresults.out','w')
summaryfile.write("Residue  Dmax   Kd \n")   

temp_Kds=[]
temp_Dmax=[]
for iy, y in enumerate(data):
    fit_params = Parameters()
    if fixedKd:
    	print "FFFFF"
    	fit_params.add( 'Kd', vary = False, value=fixedKd, min=0.0,  max=100)
    else:
    	fit_params.add( 'Kd', vary = True, value=Kd_init, min=0.0,  max=100)
    if fixDmax:
    	fit_params.add('Dmax', vary=False, value=Dmaxlist[iy])
    else:
    	fit_params.add( 'Dmax', vary=True, value=data.max(), min=0, max=data.max()*10)
#    fit_params.add( 'n_%i' % (iy+1), vary= False, value=n, min=1, max =3)
    localdata=data[iy] 
    minner = Minimizer(fcn2min, fit_params, fcn_args=(x, localdata))
    result = minner.minimize()
    final = localdata + result.residual
    report_fit(result)

    plt.plot(x_sim,simcurve(x_sim,result.params['Kd'].value,result.params['Dmax'].value))
    plt.plot(x, localdata, 'k+')
#    plt.plot(x, final, 'r')
    filename=str(outputdir)+'/fit'+str(int(residues[iy]))+'.pdf'
    titlefig="Residue "+str(int(residues[iy]))
    plt.title(titlefig)
    plt.savefig(filename)
    plt.close()
    temp_Kds.append(result.params['Kd'].value)
    temp_Dmax.append(result.params['Dmax'].value)
    simdata[iy+1][1:]=simcurve(x_sim,result.params['Kd'].value,result.params['Dmax'].value)

    summaryfile.write(str(residues[iy])+"\t "+str(result.params['Dmax'].value)+"\t "+str(result.params['Kd'].value)+"\n")   

fitted_Kds=np.array(temp_Kds)
fitted_Dmax=np.array(temp_Dmax)

summaryfile.write("Median of fitted Kd values: "+str(np.median(fitted_Kds))+"\n")

print "Median of fitted Kds ", np.median(fitted_Kds)
print "Dmax",fitted_Dmax
print "Residues",residues

Kdoutputfile = open(str(outputdir)+'/fittedKd.out','w')
Dmaxoutputfile = open(str(outputdir)+'/fittedDmax.out','w')
simcurvesfile = open(str(outputdir)+'/fittedcurves.out','w')
np.savetxt(Kdoutputfile, fitted_Kds, delimiter=',')
np.savetxt(Dmaxoutputfile, fitted_Dmax,  delimiter=',')
simcurves=np.transpose(simdata)
np.savetxt(simcurvesfile, simcurves,  delimiter=' ',newline='\n')



Kdoutputfile.close()
Dmaxoutputfile.close()
simcurvesfile.close()
summaryfile.close()
#################################################
#################################################
