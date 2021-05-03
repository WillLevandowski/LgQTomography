README: QTomo bundle

A set of seven MATLAB scripts for inversion of Lg amplitude data.
Simultaneous joint inversions for "S"ource and "R"eceiver amplification
   factors and 2D Lg attenuation (1/Q) structure are performed.
This bundle is a cleaned version of the scripts used by Levandowski, Boyd,
   McNamara, and AbdelHameid, 2021 JGR "Crustal seismic attenuation 
   structure of the central United States and Intermountain West".

QTomo_RunScript is the wrapper.
It calls the following scripts in the following order.

QTomo_SetParams sets parameters for the inversions. Geographic parameters
   should be changed as needed.

QTomo_ReadData ingests amplitude data and station/event metadata.
Ex/ "Q_data_0_75.txt" contains 0.75 Hz data and metadata.
One line per arrival
   [event,evlat,evlong,sta,stlat,stlong,dist,time,lfq,hfq,amp]
   event: Numeric ID listing origin date/time
   evlat, evlong: epicenter
   sta: Seismic station ID, alphanumeric
   stlat, stlong: station location
   dist: Epicentral distance, km.
   time: Arrival time. Unnecessary.
   lfq, hfq: Lower, upper bounds (Hz) of portion of the displacement spectrum. 
			   Nominal frequency (ex/ 0.75 Hz) = lfq/2 + hfq/2.
   amp: RMS power of displacement spectrum from lfq to hfq. units=m.
Arrivals are culled to defined longitude, latitude, and distance limits.
Each event and station must have "minv" associated arrivals or are culled.
Although not necessary, the average 1/Q and S are computed and their effects
    temporarily subtracted from the data. Inversions then run on residuals.

QTomo_traceLg is analogous to a ray tracer. 
    It sets up the inversion geometry, including a few lines that are specific 
      to the L,B,McN,AH study that should be changed for other use.
    Then, the rays are traced through a mesh of rectilinear cells. The entry 
      and exit points of each cell are determined, then the amplitude kernel
      for this ray segment are computed. (Kernels are computed radially rather
      than as annuli about each point for ease. The raypath-integrated kernels 
      are the same for annuli and radial point-kernels.)
~~~~~~~~QTomo_traceLgRays does the same thing for discrete raypaths rather than
~~~~~~~~the time-window adjoint travel time stuff used by L,B,McN,AH. 
    
QTomo_BuildSmoothingKernel builds the smoothing (penalty) function.
   The weight "smooth_coefficient" and exponent "smooth_f_exponent" were 
      determined by Tikhonov regularization of each frequency independently.
      For this dataset, the optimal weights were fit to w_optimal=40*f^0.8
      
QTomo_Inversion conducts the inversion GG*model=dd ... model ~=~ GG^-1*dd
    The model vector contains S, R, and 1/Q values. 
    The data vector is amplitude data (here, relative to constant S & 1/Q)
       and zeros attendant to smoothing matrix
    GG matrix has four components:  Binary matrix of source indices.
    								Binary matrix of receiver indices.
    								Amplitude kernels in 1/Q
    								Smoothing matrix "smooth_A"
	After one inversion, 2D Q structure is built and checked. 
	   Because smoothing/damping are ineffective as model-->0, some Q values
	      may approach infinity or even become negative.
	   These spuriously high Q are manually floored to an objective maximum.
	      Because small d(1/A) --> large d(Q) as A-->0, this step has minor
	         impact on amplitudes.
	The residuals are recomputed with the limited-Q structure.
	Additional weighting terms (like a logarithmic barrier function) are 
	    computed, such that perturbations to already-high Q are
	    hopefully damped into non-existence. Smoothing weights of 
	    average-Q and low-Q nodes are relaxed to allow tighter fit.
	Another inversion takes place to conform to the upper bounds placed on Q.
	
QTomo_WriteResults formats and outputs the results.
	Any application-specific formatting is probably for GMT.

make_contour has been a big time-saver to quickly visualize results.
	[x coordinate, y coordinate, value, pixel spacing in same units as x and y]
	It is written to be somewhat flexible.
	   Works with Lon/Lat (not to equal area), m, km, etc.
	The automated colormap editing in this script is designed specifically for
	   "value" vectors centered at or near 0 (e.g., density perturbations 
	   relative to the average at a depth, velocity pertubations, etc.), so 
	   there can be some surprises with non-zero-centered "value" vectors like Q.
	   

	
