def FitPowerlaw(fitData):
	returnVal = 0
	import pylab as plb
	import numpy as np
	import matplotlib.pyplot as plt
	from scipy.optimize import curve_fit
	from scipy import asarray as ar,exp

	n = len(fitData)
	x = 0:(n-1)
	y = fitData
	
	a = 1
	b= -1
	c = max(y)	
			
	def powerlaw(x,a,b,c):
		return a*x^b + c
		
	popt,pcov = curve_fit(powerlaw,x,y,p0=[a,b,c])

	
	plt.clf()
	plt.plot(x,y,'b+:',label='data')
	xFit = np.linspace(min(x), max(x), 5*len(x))
	plt.plot(xFit,gaussOnGradient(xFit,*popt),'r-',label='fit')
	plt.legend()
	plt.title('Peak Profile Fitting')
	plt.xlabel('R (recip. Angstroms)')
	plt.ylabel('Intensity')
	
	fileName = '..\\..\\..\\..\\peakFit_'
	meanToString = '%.4f' %  (popt[1]) 
	fileName += meanToString.replace(".", "")
	#fileName += mean
	fileName += '.png'
	print fileName
	
	figure = plt.gcf()
	figure.set_size_inches(8, 6)
	plt.savefig(fileName, bbox_inches='tight', dpi = 150)
	plt.show(block=False)
	
	#Need to concatenate popt numpy array into string??
	return popt[2]
	
