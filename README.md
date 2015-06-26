#BALparams.py
Author: Jesse A. Rogerson - jesserogerson.com

The code implementation was written by the above, but with helpful (and crucial) contributions from: **Catherine J. Grier, Patrick B. Hall, Daniel E. Vandenberk**

-----

### Synopsis

This code calculates the 'BALnicity Index' of a quasar given some user-supplied spectrum.

In the literature, the definition of the BALnicity Index has gone through many iterations. Here we use BALnicity Index as a catch-all term for any measurement in the same spirit as the original BI (see below). This code is meant to be a one-stop-shop for all definitons. It allows the user to either:

a. choose from a set of predefined indexes  
b. manually decide values of some or all of the parameters to define their own BALnicity Index

###Predefined Indexes Available:

1. BI - Weymann et al. 1991, ApJ, 373, 23
2. BI_0 - Gibson et al., 2008, ApJ, 675, 985
3. AIT - Trump et al., 2006, ApJS, 165, 1
4. AI450 - Hall et al. 2002, ApJS, 141, 267
5. DI - Paris et al. 2012, A&A, 548, 66


###For Help:
`$> ./BALparams.py -h`

This will output the various REQUIRED or 'positional arguments' as well as the various 'optional arguments.'

###Example1 (The Default Execution):
The only required inputs (or positional arguments) from the user are:

**file**: An ASCII spectrum of a quasar. It should be normalized to avoid bad pixels. It should also be in the rest-frame. The first three columns must be: #wavelength flux flux_err

**zem**:  The Redshift of the quasar

With these, the user can execute:

`$> ./BALparams.py file zem`

Upon executing without any other parameters specified,
BALparams.py will default to the *original* BI from
Weymann et al. 1991, ApJ, 373, 23

###Example2 (Using Predefined Indexes):

There are a few pre-set BALnicity indexes that are available in this script: BI, BI_0, AI450, DI

If the user specifies one of these via the '-index' parameter on the commandline, it will automatically set the code to calculate that index.

`$> ./BALparams.py file zem -index AI450`

In the above example, all values of velocity ranges, minimum velocity, and more will be set to those used in the paper where AI450 was defined.

###Example3 (Fully Manual):

The user may choose to define any (or all) of the 'optional arguments' in the 'help' section. However, in order to make this work, the user MUST specify '-index man,' (see below)

e.g.,

`$> ./BALparams.py file zem -index man -vlo 0 -vhi 35000`

`$> ./BALparams.py file zem -index man -vhi 35000 -zerr 0.056 -v 1000 -inc y`

And any 'optional parameters' you do not set will remain the default BI values from Weymann et al. 1991

###Notes:

a. If you set '-index' to be one of the available pre-set indexes, it will override any additional optional parameters you set.

b. The code does not validate your optional parameters. So if you define a '-vlo' that is higher than the defined '-vhi,' the code will NOT yell at you.

c. The code will print out the configuration of the
variables each time it runs, so you will know what you set.

### Example Plot:

Along with the script itself, an example output plot is provided named 'BALplot.eps' and included here. The plot was created using the following command:

`
$> BALparams.py spectrum.dat 2.5 -index man -vlo 0 -vhi 60000 -v 450 -inc y -out BALplot.png`

![BALplot](BALplot.png =600x)
The shaded in blue regions are the portions of the spectrum that meet the BALnicity criteria as stipulated on the command line.


### Installation

No Installation needed, just download and execute script.

### License

see LICENSE.txt
