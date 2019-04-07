import ellc
import emcee
import pyfits
import numpy as np
import scipy.optimize as opt
import corner
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator

#Returns an array containing the quaterly light curves of the planet
def open_files(planet):
	#KOI-1741
	if str(planet) == 'KOI-1741':
		Q4 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2010078095331_llc.fits')
		Q5 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2010174085026_llc.fits')
		Q6 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2010265121752_llc.fits')
		Q7 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2010355172524_llc.fits')
		Q8 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2011073133259_llc.fits')
		Q9 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2011177032512_llc.fits')
		Q10 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2011271113734_llc.fits')
		Q11 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2012004120508_llc.fits')
		Q12 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2012088054726_llc.fits')
		Q13 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2012179063303_llc.fits')
		Q14 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2012277125453_llc.fits')
		Q15 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2013011073258_llc.fits')
		Q16 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2013098041711_llc.fits')
		Q17 = pyfits.open('/home/astro/phunbw/Documents/James_URSS/Light curves/kplr008180020-2013131215648_llc.fits')
		
	Q = [Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17]
	return Q

#Extracts the data, cleans and normalises - Returns Time, Flux, Error arrays
def obtain_arrays(planet):
	#Obtain the quarters of data
	Q = open_files(str(planet))
	
	#Define Time, Flux and Error arrays
	Time = []
	Flux = []
	Error = []
	
	for i in Q:
		#Extract data
		dat = i[1].data
		x_i = dat['Time']
		y_i = dat['PDCSAP_FLUX']
		quality_i = dat['SAP_QUALITY']
		error_i = dat['PDCSAP_FLUX_ERR']
		i.close()
    
		#Clean data
		Boolean = quality_i == 0  
		x_i = x_i[Boolean]
		y_i = y_i[Boolean]
		error_i = error_i[Boolean]
    
		#Remove nans from x_i
		nan_check = ~np.isnan(x_i)
		x_i = x_i[nan_check]
		y_i = y_i[nan_check]
		error_i = error_i[nan_check]
    
		#Remove nans from y_i
		nan_check = ~np.isnan(y_i)
		x_i = x_i[nan_check]
		y_i = y_i[nan_check]
		error_i = error_i[nan_check]
		
		#Remove nans from error_i
		nan_check = ~np.isnan(error_i)
		x_i = x_i[nan_check]
		y_i = y_i[nan_check]
		error_i = error_i[nan_check]
    
		#Obtain median for relative flux - normalise
		median_i = np.median(y_i)
		y_i /= median_i
		error_i /= median_i
    
		#Obtain arrays
		Time = np.append(Time, x_i) 
		Flux = np.append(Flux, y_i)
		Error = np.append(Error, error_i)
		
	return Time, Flux, Error

#Removes the eclipses from Time, Flux and Error ready for polyfit
def remove_eclipses(planet, t0, period):	
	Time, Flux, Error = obtain_arrays(planet)
    
	number = 0
	count = t0

	#Remove the primary and secondary eclipses
	while count < max(Time):
        
		count = t0 + number*period
    
		#Define the start and end of the relevant primary eclipse
		pri_eclipse_start = count - 0.18
		pri_eclipse_end = count + 0.18
    
		#Remove the primary eclipses
		Boolean = (Time < pri_eclipse_start) | (Time > pri_eclipse_end)

		Time = Time[Boolean]
		Flux = Flux[Boolean]
		Error = Error[Boolean]
    
		#Define the start and end of the relevant secondary eclipse
		sec_eclipse_start = count + 0.48736*period
		sec_eclipse_end = count + 0.52492*period
    
		#Remove the secondary eclipses
		Boolean = (Time < sec_eclipse_start) | (Time > sec_eclipse_end)
    
		Time = Time[Boolean]
		Flux = Flux[Boolean]
		Error = Error[Boolean]
    
		number += 1
    
	return Time, Flux, Error

#Fits polynomials to chunks of the remaining dataset for flattening 
def polynomial_fit(planet, t0, period, step_size):  
    Time, Flux, Error = obtain_arrays(planet)
    
    #Store original Time, Flux and Error arrays for future use
    Time_with_eclipses = Time
    Flux_with_eclipses = Flux 
    Error_with_eclipses = Error   
    
    Time, Flux, Error = remove_eclipses(planet, t0, period)
    
    epoch = min(Time)
    number = 0
    count = epoch
    step_size = float(step_size)

    predicted_final = []
    Time_flat = []
    Flux_for_plot = []
    Error_for_plot = []

    while count < (max(Time) + step_size):
		#Define start and end of chunk
        start = epoch + number*step_size
        end = epoch + number*step_size + step_size
		
		#Remove all the data outside of the chunk
        Boolean = (Time > start) & (Time < end)

        Time_chunk = Time[Boolean]
        Flux_chunk = Flux[Boolean]
        Error_chunk = Error[Boolean]
		
		#Repeat with the original array to include the eclipses
        Boolean = (Time_with_eclipses > start) & (Time_with_eclipses < end)
    
        Time_chunk_with_eclipses = Time_with_eclipses[Boolean]
        Flux_chunk_with_eclipses = Flux_with_eclipses[Boolean]
        Error_chunk_with_eclipses = Error_with_eclipses[Boolean]
		
		#Fit polynomial to chunk, divide data by model to flatten
        if len(Time_chunk) > 0:
            model = np.polyfit(Time_chunk, Flux_chunk, 3)
            predicted = np.polyval(model, Time_chunk_with_eclipses)
        
            Time_flat = np.append(Time_flat, Time_chunk_with_eclipses)
            Flux_for_plot = np.append(Flux_for_plot, Flux_chunk_with_eclipses)
            Error_for_plot = np.append(Error_for_plot, Error_chunk_with_eclipses)
            predicted_final = np.append(predicted_final, predicted)
    
        number += 1
        count = epoch + number*step_size
    
    Flux_flat = Flux_for_plot / predicted_final
    Error_flat = Error_for_plot / predicted_final
    
    return Time_flat, Flux_flat, Error_flat

#For plotting the flattened light curve    
def plot_flatten(planet, t0, period, step_size):  
    Time_flat, Flux_flat, Error_flat = polynomial_fit(planet, t0, period, step_size)
    Time, Flux, Error = obtain_arrays(planet)
    
    plt.plot(Time_flat, Flux_flat, 'b')
    plt.plot(Time, Flux, 'r')

    plt.title(str(planet) + ' Flattened Light Curve')
    plt.xlabel('Time [BJD-2454833]')
    plt.ylabel('Relative Flux')

    plt.show()

#Prepares the necessary arrays for plotting the phase curve of the planet  
def phase_curve(planet, t0, period, step_size):	
	Time_flat, Flux_flat, Error_flat = polynomial_fit(planet, t0, period, step_size)
	
	#Find phase
	time_difference = Time_flat - t0
	Phase = np.mod(time_difference, period) / period

	#Sort the arrays into order
	Phase_index = np.argsort(Phase)	
	Phase = Phase[Phase_index]
	Flux_flat = Flux_flat[Phase_index]
	Error_flat = Error_flat[Phase_index]
	
	#Histogram and digitize to get bins, then calculate mean of fluxes in each bin
	bin_error = []
	bins = np.histogram(Phase, 500)[1]
	n = np.histogram(Phase, 500)[0]
	binned = np.digitize(Phase, bins)
	
	bin_phases = np.array([Phase[binned == i].mean() for i in range(1, len(bins))])
	bin_means = np.array([Flux_flat[binned == i].mean() for i in range(1, len(bins))])
	bin_means_err = np.array([Error_flat[binned == i].mean() for i in range(1, len(bins))])
	
	for i in range(len(bin_means_err)):
		error = bin_means_err[i] / np.sqrt(n[i])
		bin_error = np.append(bin_error, error)
	
	return Phase, Phase_index, bin_phases, bin_means, bin_error

#Removes all the data outside of the eclipses to speed up ellc fitting   
def keep_eclipses(planet, t0, period, step_size):	
	Time_flat, Flux_flat, Error_flat = polynomial_fit(planet, t0, period, step_size)
	
	number = 0
	count = t0
	
	Time_final = []
	Flux_final = []
	Error_final = []

	#Remove data outside eclipses
	while count < max(Time_flat):
        
		count = t0 + number*period
    
		#Define the start and end points of the desired chunks
		pri_eclipse_start = count - 0.05*period
		pri_eclipse_end = count + 0.05*period
		sec_eclipse_start = count + 0.45*period
		sec_eclipse_end = count + 0.55*period
    
		#Preserve the primary eclipses
		Boolean = ((Time_flat > pri_eclipse_start) & (Time_flat < pri_eclipse_end)) | ((Time_flat > sec_eclipse_start) & (Time_flat < sec_eclipse_end))
    
		Time_eclipses = Time_flat[Boolean]
		Flux_eclipses = Flux_flat[Boolean]
		Error_eclipses = Error_flat[Boolean]
		
		if len(Time_eclipses) > 0:
			Time_final = np.append(Time_final, Time_eclipses)
			Flux_final = np.append(Flux_final, Flux_eclipses)
			Error_final = np.append(Error_final, Error_eclipses)
    
		number += 1
    
	return Time_final, Flux_final, Error_final

#Takes in the binary parameters and returns an ellc fit for the light curve	
def light_curve_model(t_obs, t0, period, radius_1, radius_2, sbratio, incl, f_c, f_s, ldc_1):
	
	lc_model = ellc.lc(t_obs=t_obs,radius_1=radius_1,radius_2=radius_2,sbratio=sbratio,incl=incl,t_zero=t0,period=period,q=0.21864444,a=14.6,t_exp=0.02044,n_int=5,f_c=f_c,f_s=f_s,ldc_1=ldc_1,ldc_2=0.6,ld_1='lin',ld_2='lin',grid_1='very_sparse',grid_2='very_sparse')
	return lc_model

#Takes in the binary parameters and returns an ellc fit for the radial velocity curve
def rv_curve_model(t_obs, t0, period, f_c, f_s, a, offset, ):
	
	rv1,rv2 = ellc.rv(t_obs=t_obs,radius_1=in_radius_1,radius_2=in_radius_2,sbratio=in_sbratio,incl=in_incl,t_zero=t0,period=period,q=0.21864444,a=a,f_c=f_c,f_s=f_s,grid_1='very_sparse',grid_2='very_sparse')
	
	#Account for instrumental offset
	rv1 = rv1 - offset
	return rv1
	
#Initial guesses of the parameters
in_radius_1 = 1.429175874 / 16     #solar radii
in_radius_2 = 0.3204153338 / 16    #solar radii
in_sbratio = 0.0311123506714
in_incl = 89.6073820235
in_t0 = 356.620105511              #from light curve
in_period = 5.80312668749
#in_q = 0.21864444
#in_a =                        #solar radii
in_f_c = 0.107117375335
in_f_s = -0.019639617129
in_ldc_1 = 0.459650163699
#in_ldc_2 = 
in_offset = 56.1958919779

initial = [in_radius_1, in_radius_2, in_sbratio, in_incl, in_t0, in_period, in_f_c, in_f_s, in_a, in_offset, in_ldc_1]

############
### MCMC ###
############

#Extract the necessary data	
Time_final, Flux_final, Error_final = keep_eclipses(planet='KOI-1741', t0=in_t0, period=in_period, step_size=2)

#Define the relevant arrays

x_rv = [1443.82210,1446.92356,1457.82103,1459.79856,1460.90038,1461.91270]
x_rv = np.array(x_rv)
x_rv = x_rv + 167.0
y_rv = [-74.9556,-33.3380,-51.6407,-38.7746,-65.6135,-81.6842]
y_rv = np.array(y_rv)
yerr_rv = [1.180,0.341,0.167,0.182,0.122,0.613]
yerr_rv = np.array(yerr_rv)

x_lc = Time_final
y_lc = Flux_final
yerr_lc = Error_final

def lnprior(theta):
	radius_1, radius_2, sbratio, incl, t0, period, f_c, f_s, a, offset, ldc_1 = theta
	#Uniform priors for the parameters in theta
	if 0.5/14.6 < radius_1 < 2.0/14.6 and 0.1/14.6 < radius_2 < 0.5/14.6 and 0.01 < sbratio < 0.5 and 70 < incl < 100 and 356.600 < t0 < 356.650 and 5.79 < period < 5.81 and -0.8 <= f_c < 0.8 and -0.8 <= f_s < 0.8 and 14.0 < a < 18.0 and 0.0 < offset < 100.0 and 0.3 < ldc_1 < 0.7:
		return 0.0
	return -np.inf
	
def lnlike(theta, x_lc, y_lc, yerr_lc, x_rv, y_rv, yerr_rv):
	radius_1, radius_2, sbratio, incl, t0, period, f_c, f_s, a, offset, ldc_1 = theta
	
	#Light curve likelihood function
	model_lc = light_curve_model(t_obs=x_lc, t0=t0, period=period, radius_1=radius_1, radius_2=radius_2, sbratio=sbratio, incl=incl, f_c=f_c, f_s=f_s, ldc_1=ldc_1)
	if (True in np.isnan(model_lc)) or np.min(model_lc) <= 0:
		lnlike_lc = -np.inf
	else:
		inv_sigma2_lc = 1.0/(yerr_lc**2)
		lnlike_lc = -0.5*(np.sum((y_lc-model_lc)**2*inv_sigma2_lc - np.log(inv_sigma2_lc)) - np.log(len(y_lc) + 1))
	
	#Rv curve likelihood function
	model_rv = rv_curve_model(t_obs=x_rv, t0=t0, period=period, a=a, f_c=f_c, f_s=f_s, offset=offset)
	if (True in np.isnan(model_rv)):
		lnlike_rv = -np.inf
	else:
		inv_sigma2_rv = 1.0/(yerr_rv**2)
		lnlike_rv = -0.5*(np.sum((y_rv-model_rv)**2*inv_sigma2_rv - np.log(inv_sigma2_rv)) - np.log(len(y_rv) + 1))
	
	#Sum to get overall likelihood function
	lnlike = lnlike_lc + lnlike_rv
	return lnlike
	
def lnprob(theta, x_lc, y_lc, yerr_lc, x_rv, y_rv, yerr_rv):
	lp = lnprior(theta)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnlike(theta, x_lc, y_lc, yerr_lc, x_rv, y_rv, yerr_rv)
	
#Set up the sampler
ndim, nwalkers, nsteps = 11, 50, 5000
weights = [5e-4, 5e-4, 5e-3, 1e-2, 1e-4, 5e-5, 5e-3, 5e-3, 1e-2, 1e-2, 1e-2]
pos = [initial + weights*np.random.randn(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x_lc, y_lc, yerr_lc, x_rv, y_rv, yerr_rv))

#Run the production chain
print("Running MCMC...")
sampler.run_mcmc(pos, nsteps, rstate0=np.random.get_state())
print("Done.")

########################################################
### Plot and save the times series of each parameter ###
########################################################

plt.plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
plt.axhline(in_radius_1, color="#888888", lw=2)
plt.ylabel('radius_1')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_radius_1.png')
plt.clf()

plt.plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
plt.axhline(in_radius_2, color="#888888", lw=2)
plt.ylabel('radius_2')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_radius_2.png')
plt.clf()

plt.plot(sampler.chain[:, :, 2].T, color="k", alpha=0.4)
plt.axhline(in_sbratio, color="#888888", lw=2)
plt.ylabel('sbratio')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_sbratio.png')
plt.clf()

plt.plot(sampler.chain[:, :, 3].T, color="k", alpha=0.4)
plt.axhline(in_incl, color="#888888", lw=2)
plt.ylabel('incl')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_incl.png')
plt.clf()

plt.plot(sampler.chain[:, :, 4].T, color="k", alpha=0.4)
plt.axhline(in_t0, color="#888888", lw=2)
plt.ylabel('t0')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_t0.png')
plt.clf()

plt.plot(sampler.chain[:, :, 5].T, color="k", alpha=0.4)
plt.axhline(in_period, color="#888888", lw=2)
plt.ylabel('period')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_period.png')
plt.clf()

plt.plot(sampler.chain[:, :, 6].T, color="k", alpha=0.4)
plt.axhline(in_f_c, color="#888888", lw=2)
plt.ylabel('f_c')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_f_c.png')
plt.clf()

plt.plot(sampler.chain[:, :, 7].T, color="k", alpha=0.4)
plt.axhline(in_f_s, color="#888888", lw=2)
plt.ylabel('f_s')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_f_s.png')
plt.clf()

plt.plot(sampler.chain[:, :, 8].T, color="k", alpha=0.4)
plt.axhline(in_a, color="#888888", lw=2)
plt.ylabel('a')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_a.png')
plt.clf()

plt.plot(sampler.chain[:, :, 9].T, color="k", alpha=0.4)
plt.axhline(in_offset, color="#888888", lw=2)
plt.ylabel('offset')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_offset.png')
plt.clf()

plt.plot(sampler.chain[:, :, 10].T, color="k", alpha=0.4)
plt.axhline(in_ldc_1, color="#888888", lw=2)
plt.ylabel('ldc_1')
plt.xlabel('step number')
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_ldc_1.png')
plt.clf()

###################################################
### Calculate the most likely set of parameters ###
###################################################

#Get user to input the burnin period, after they have seen the time series of each parameter
burnin = int(raw_input('Enter burnin period: '))
samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))

#Most likely set of parameters
radius_1 = np.median(samples[:,0])
radius_2 = np.median(samples[:,1])
sbratio = np.median(samples[:,2])
incl= np.median(samples[:,3])
t0 = np.median(samples[:,4])
period = np.median(samples[:,5])
f_c = np.median(samples[:,6])
f_s = np.median(samples[:,7])
a = np.median(samples[:,8])
offset = np.median(samples[:,9])
ldc_1 = np.median(samples[:,10])

print 'radius_1 = ',radius_1 ,' +/- ',np.std(samples[:,0])
print 'radius_2 = ',radius_2 ,' +/- ',np.std(samples[:,1])
print 'sbratio = ',sbratio ,' +/- ',np.std(samples[:,2])
print 'incl = ',incl ,' +/- ',np.std(samples[:,3])
print 't0 = ',t0 ,' +/- ',np.std(samples[:,4])
print 'period = ',period ,' +/- ',np.std(samples[:,5])
print 'f_c = ',f_c ,' +/- ',np.std(samples[:,6])
print 'f_s = ',f_s ,' +/- ',np.std(samples[:,7])
print 'a = ',a ,' +/- ',np.std(samples[:,8])
print 'offset = ',offset ,' +/- ',np.std(samples[:,9])
print 'ldc_1 = ',ldc_1 ,' +/- ',np.std(samples[:,10])

##########################
### Plot triangle plot ###
##########################

fig = corner.corner(samples, labels=["$radius_1$", "$radius_2$", "$sbratio$", "$incl$", "$t0$", "$period$", "$f_c$", "$f_s$", "$a$", "$offset$", "$ldc_1$"], truths=initial, plot_contours=False)
plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_corner.png')
fig.clf()

##############################################################
### Take most likely set of parameters and plot the models ###
##############################################################

#Prep for plotting
t_model = np.linspace(1605,1635,10000)

#Calculate final models
lc_model = light_curve_model(x_lc, t0, period, radius_1, radius_2, sbratio, incl, f_c, f_s, ldc_1)
rv_model = rv_curve_model(t_model, t0, period, f_c, f_s, a, offset)

#Phase the original flattened data using the calculated t0 and period values
Phase, Phase_index, bin_phases, bin_means, bin_error = phase_curve(planet='KOI-1741', t0=t0, period=period, step_size=2)

#Subplot 1 - For the light curve
fig = plt.figure()
gs = gridspec.GridSpec(3,3)

ax1 = fig.add_subplot(gs[0,0:])
plt.errorbar(bin_phases, bin_means, yerr=bin_error, fmt='r.')
plt.plot(Phase, lc_model[Phase_index], 'b')
plt.errorbar(bin_phases-1, bin_means, yerr=bin_error, fmt='r.')
plt.plot(Phase-1, lc_model[Phase_index], 'b')

plt.xlabel('Phase')
plt.ylabel('Relative Flux')

#Subplot 2 - For the primary eclipse (zoom in)
ax2 = fig.add_subplot(gs[1,0])
plt.errorbar(bin_phases, bin_means, yerr=bin_error, fmt='r.')
plt.plot(Phase, lc_model[Phase_index], 'b')
plt.errorbar(bin_phases-1, bin_means, yerr=bin_error, fmt='r.')
plt.plot(Phase-1, lc_model[Phase_index], 'b')

plt.xlabel('Phase')
plt.ylabel('Relative Flux')

#Subplot 3 - For the secondary eclipse (zoom in)
ax3 = fig.add_subplot(gs[2,0])
plt.errorbar(bin_phases, bin_means, yerr=bin_error, fmt='r.')
plt.plot(Phase, lc_model[Phase_index], 'b')
plt.errorbar(bin_phases-1, bin_means, yerr=bin_error, fmt='r.')
plt.plot(Phase-1, lc_model[Phase_index], 'b')

plt.xlabel('Phase')
plt.ylabel('Relative Flux')

#Subplot 4 - For the radial velocities
ax4 = fig.add_subplot(gs[1:,1:])
plt.plot(t_model, rv1,'b')
plt.errorbar(x_rv, y_rv, yerr=yerr_rv, fmt='r.')

plt.xlabel('Time [BJD-2454833]')
plt.ylabel('rv1 [km/s]')

plt.savefig('/home/astro/phunbw/Documents/James_URSS/Figures/emcee both '+str(nsteps)+' steps/emcee_both_'+str(nsteps)+'steps_'+str(nwalkers)+'walkers_models.png')
fig.clf()

