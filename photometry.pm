package photometry;
my $MODVERSION = '$Id: photometry.pm,v 1.12 2008/11/12 23:27:09 sgroom Exp $';


##########################################################################
# NASA Exoplanet Science Institute (NExScI) -- http://nexsci.caltech.edu #
#       Copyright (c) 1997-2008 California Institute of Technology       #
#                     U.S. Government Sponsorship                        #
# Under contract with the National Aeronautics and Space Administration  #
##########################################################################


## photometry.pm -- define (hopefully) useful photometric quantities
## like wavelengths, passbands, zeropoints, and the like

use strict;
use warnings;

require 5.000;

use Exporter();
our @ISA = qw(Exporter);
our @EXPORT = qw(fLambda fNu bandParameters matchBand fLambdaToFNu fNuToFLambda formatPhotRec
	UNITS_MAG UNITS_FLAMBDA UNITS_FNU);

use miscString;
use gcConfig;

my $cfg = gcConfig::config();

my $center = 0;
my $passband = 1;
my $zeroFLambda = 2;
my $zeroFNu = 3;

use constant UNITS_MAG => 0;		# Report in units of Mag
use constant UNITS_FLAMBDA => 1; 	# Report in units of fLambda (erg cm^-2 s^-1 um^-1)
use constant UNITS_FNU => 2;		# Report in units of fNu (Jy;  10^-26 W m^-2 Hz^-1)

### Johnson system hash
my %Johnson = (
	    "U" => [360,68,4.22e-5,1.823e+3],		## Vis photometry properties from AQ table 15.6 (p 387)
	    "B" => [450,98,6.40e-5,4.130e+3],		## Flambda zeropoints in erg cm^-2 s^-1 um^-1
	    "V" => [555,89,3.75e-5,3.781e+3],		## Fnu units in Jy (10^-26 W m^-2 Hz^-1)
	    "R" => [693,210,1.87e-5,2.89e+3],		## New zeropts for R and I from
##	    "R" => [670,220,1.74e-5,2.941e+3],		##  Fukugita et al 1995 table 9
	    "I" => [879,171,9.12e-6,2.28e+3],		# also cf. http://www.astro.utoronto.ca/~patton/astro/mags.html#conversions
##	    "I" => [900,240,8.4e-6,2.276e+3],

	    "J" => [1215,260,3.31e-6,1630],     ## IR photometry properties from AQ table 7.5 (p 150)
	    "H" => [1654,290,1.15e-6,1050],     ## Flambda zeropoints in erg cm^-2 s^-1 um^-1
	    "K" => [2179,410,4.14e-7,655],      ## Fnu units in Jy (10^-26 W m^-2 Hz^-1)
	    "L" => [3500,700,6.59e-8,276],
	    "M" => [4769,450,2.11e-8,160],
	    "N" => [10472,5190,9.63e-10,35.2],
	    "Q" => [20130,7800,7.18e-11,9.70],
	    );

my %KronComet = (
     "OH" => [309.7,5.8,0.00010560,3377],
     "NH" => [336.1,5.4,0.00008420,3171],
     "UVc" => [344.9,7.9,0.00007802,3094],
     "CN" => [386.9,5.6,0.00008600,4292],
     "C3" => [406.3,5.8,0.00008160,4491],
     "COp" => [426.6,6.4,0.00007323,4443],
     "Bc" => [445.3,6.1,0.00006210,4105],
     "C2" => [513.5,11.9,0.00003887,3417],
     "Gc" => [525.9,5.6,0.00003616,3334],
     "H2Op" => [702.8,16.4,0.00001380,2272],
     "Rc" => [713.3,5.8,0.00001316,2232],
     );

my %Oja = (
     "m45" => [450.8,5.0,0.00006221,4215],
     "m42" => [426.9,4.3,0.00013193,8016],
     "m41" => [417.6,4.0,0.00013867,8062],
     );

	 
my %Eggen102 = (								## interpolated using "system updating 120120a.xls" spreadsheet
#     "F62" => [624.0,27.6,0.00002411,3127],		## http://obswww.unige.ch/gcpd/ph15.html
#     "F65" => [650.6,27.2,0.00002128,2986],
#     "F102" => [1018.4,32.0,0.00000558,2244],
     "F62" => [624.0,27.6,0.00002507,3254],
     "F65" => [650.6,27.2,0.00002166,3056],
     "F102" => [1018.4,32.0,0.00000450,1556],
     );

### Cousins
my %Cousins = (
	    "Rc" => [659,157,2.15e-5,3.020e+3],  ## From Fukugita et al 1995 table 9
	    "Ic" => [806,154,1.10e-5,2.380e+3],
	    );


### Geneva photometry
## Best reference for Geneva photometric system is
## Rufener and Nicolet 1988, A&A 206, 357-374
my %Geneva = (
	   "U" =>  [346.4,15.9,5.7544e-5,2304.092976],
	   "B1" => [401.5,18.8,6.76083e-5,3641.665066],
	   "B" =>  [422.7,28.2,2.88403e-5,1717.51287],
	   "B2" => [447.6,16.3,0.000108443,7259.387846],
	   "V1" => [539.5,20.2,7.46105e-5,7252.704781],
	   "V" =>  [548.8,29.6,3.73594e-5,3754.049545],
	   "G" =>  [580.7,20.0,9.48855e-5,10697.93611],
	   );


### Stromgren (uvby) photometry
## zeropoints from Gray 1998, AJ 116, 482.
my %Stromgren = (
	      "u" => [350,34,11.72e-5,4.79e+03],
	      "v" => [411,20,8.66e-5,4.88e+03],
	      "b" => [467,16,5.89e-5,4.28e+03],
	      "y" => [547,24,3.73e-5,3.72e+03],
	      );


## 2Mass photometry
## Note that 2Mass bands are J,H, and Ks
## zeropoints from Cohen, Wheaton, and Megeath 2003 AJ 126, 1090
## (see also http://ssc.spitzer.caltech.edu/tools/magtojy/ref.html)
## Updated 30 Jan 2004 AB
my %TwoMass = (
	    "J" => [1235,162,3.14e-6,1594],      ## Cohen, Wheaton, and Megeath 2003
	    "H" => [1662,251,1.11e-6,1024],      ## Flambda zeropoints in erg cm^-2 s^-1 um^-1
	    "Ks" => [2159,262,4.29e-7,666.8],    ## Fnu units in Jy (10^-26 W m^-2 Hz^-1)
	    );


## CIT system photometry from Elias et al 1982 AJ 87, 1029
## CIT IR system bands are J, H, K, and L:
## J = 1.13 -- 1.37 um
## H = 1.50 -- 1.80 um
## K = 2.01 -- 2.42 um
## L = 3.22 -- 3.76 um

my %CIT = (
	"J" => [1250,240,3.14e-6,1594],  ## These zeropoints are probably bogus
	"H" => [1625,300,1.11e-6,1024],  ## looking for a good reference!!!
	"K" => [2217,410,4.29e-7,666.8],
	"L" => [3490,540,6.59e-8,276],
	);


## Walraven system photometry from Gerard van Belle
## W = 354.3 nm/48 nm
## B = 435.5 nm/88.5 nm
## V = 543.5 nm/86 nm
## R = 694.5 nm/116 nm

my %Walraven = (
#	"W" => [354.3,48,3.42e-5,1433],    ## These zeropoints are from GvB
#	"B" => [435.5,88.5,6.08e-5,3846],
#	"V" => [543.5,86,3.76e-5,3701],
#	"R" => [694.5,116,1.82e-5,2920],
#     "W" => [354.5,48.0,0.00003227,1352],	# updated zeropoints 2012-08-01
#     "B" => [435.5,88.5,0.00005986,3785],	# based on Glushneva 1998b
#     "V" => [543.5,86.0,0.00003572,3518],
#     "R" => [694.5,116.0,0.00001715,2758],
     "W" => [354.5,48.0,0.00003272,1371],	# updated zeropoints 2012-08-02
     "B" => [435.5,88.5,0.00006024,3809],	# based on Glushneva 1998b and Burnashev 1985
     "V" => [543.5,86.0,0.00003616,3561],
     "R" => [694.5,116.0,0.00001745,2807],
	);


## Hipparcos and Tycho photometry
my %Hipparcos = (
	     "Hp" => [500,200,3,4],  ## Real values not yet defined
	     );

my %Tycho = (
	 "Bt" => [428,50,6.40e-5,4.130e+3],  ## Taken from Scales et al 1992 and Grobmann et al 1995
	 "Vt" => [534,100,3.75e-5,3.781e+3],  ## Real zeropoint values not yet found!!!
	 );


## Thuan-Gunn system - from Fukujita Table 9
my %ThuanGunn = (
	 "u" => [353.6,41.2,3.33e-5,1.38e+3],
	 "v" => [399.2,46.9,6.62e-5,3.50e+3],
	 "g" => [492.7,70.9,4.84e-5,3.89e+3],
	 "r" => [653.8,89.3,2.09e-5,2.96e+3],
	 );

## SDSS system - from Fukujita Table 9
my %SDSS = (
	 "u" => [358.5, 55.6,3.67e-5, 1.54e+3],
	 "g" => [485.8,129.7,5.11e-5, 3.93e+3],
	 "r" => [629.0,135.8,2.40e-5, 3.12e+3],
	 "i" => [770.6,154.7,1.28e-5, 2.51e+3],
	 "z" => [922.2,153.0,0.783e-5,2.19e+3],
	 );

## DDO (David Dunlap Observatory) photometric system
## See McClure 1976, AJ 81, 182 for passband definitions
## relevant Vizier docs at http://vizier.u-strasbg.fr/viz-bin/Cat?II/80
## Zeropoints are self-derived using Vega spectrum
## NOT YET SURE about zeropoints!!!
my %DDO = (
#	"F35" => [346,38.3,3.366e-05,1.35e3],
#	"F38" => [381.5,33,5.187e-05,2.53e3],
#	"F41" => [416.6,8.3,7.502e-05,4.36e3],
#	"F42" => [425.7,7.3,7.272e-05,4.41e3],
#	"F45" => [451.7,7.6,6.288e-05,4.29e3],
#	"F48" => [488.6,18.6,4.587e-05,3.66e3],
#	"F51" => [513,13.4,4.367e-05,3.84e3],
#     "m35" => [346.0,38.3,0.00011512,4595],		# from burnashev 1985 spectrophotometry
#     "m38" => [381.5,33.0,0.00005304,2574],		# derived on 2012-08-06
#     "m41" => [416.6,8.3,0.00019797,11455],		# [ ] verified?
#     "m42" => [425.7,7.3,0.00018474,11161],
#     "m45" => [451.7,7.6,0.00011662,7933],
#     "m48" => [488.6,18.6,0.00004380,3486],
     "m35" => [346.0,38.3,0.00011734,4683],
     "m38" => [381.5,33.0,0.00006448,3129],
     "m41" => [416.6,8.3,0.00020003,11574],
     "m42" => [425.7,7.3,0.00018820,11370],
     "m45" => [451.7,7.6,0.00013562,9225],
     "m48" => [488.6,18.6,0.00004999,3979],
	);

my %thirteencolor = (
#     "m33" => [337.1,9.0,0.00003348,1268],			#Calibration good for A&G Main sequence ...
#     "m35" => [353.6,9.1,0.00003205,1336],
#     "m37" => [375.1,9.3,0.00004361,2046],
#     "m40" => [403.0,20.0,0.00007064,3825],
#     "m45" => [457.1,23.9,0.00006291,4382],
#     "m52" => [518.3,23.0,0.00004331,3879],
#     "m58" => [582.7,19.7,0.00003087,3495],
#     "m63" => [635.6,26.8,0.00002263,3048],
#     "m72" => [724.1,61.7,0.00001658,2898],
#     "m80" => [800.0,46.9,0.00001209,2579],
#     "m86" => [858.4,52.0,0.00001001,2458],
#     "m99" => [983.1,59.0,0.00000708,2280],
#     "m110" => [1108.4,69.6,0.00000471,1927],
#    												# Optimizing calibration for K&M giants, trimmed up with late G,K main sequence
     "m33" => [337.1,9.0,0.00003269,1239],
     "m35" => [353.6,9.1,0.00003092,1289],
     "m37" => [375.1,9.3,0.00004375,2052],
     "m40" => [403.0,20.0,0.00006791,3677],
     "m45" => [457.1,23.9,0.00006126,4267],
     "m52" => [518.3,23.0,0.00004278,3831],
     "m58" => [582.7,19.7,0.00003072,3478],
     "m63" => [635.6,26.8,0.00002270,3057],
     "m72" => [724.1,61.7,0.00001636,2860],
     "m80" => [800.0,46.9,0.00001202,2564],
     "m86" => [858.4,52.0,0.00000992,2437],
     "m99" => [983.1,59.0,0.00000716,2307],
     "m110" => [1108.4,69.6,0.00000493,2017],
     );

## Vilnius Photometry (see Straizys+ 1989)
## Filters for Vilnius photometry are given in Vizier documentation
## for Straizys catalog http://vizier.u-strasbg.fr/viz-bin/Cat?II/157A
## Zeropoints are self-derived using Vega spectrum
## NOTE that Vilnius mags for Vega are not zero (not even close!!! U=1.93!!!)
my %Vilnius = (
	    "U" => [345,40,2.00e-04,7.96e3],
	    "P" => [374,26,1.42e-04,6.65e3],
	    "X" => [405,22,1.07e-04,5.86e+3],
	    "Y" => [466,26,6.73e-05,4.89e+3],
	    "Z" => [516,21,4.61e-05,4.10e+3],
	    "V" => [544,26,3.78e-5,3.74e+3],
	    "S" => [655,22,1.81e-05,2.60e+3],
	    );

my %Straizys = (							# same as Vilnius
	    "U" => [345,40,2.00e-04,7.96e3],
	    "P" => [374,26,1.42e-04,6.65e3],
	    "X" => [405,22,1.07e-04,5.86e+3],
	    "Y" => [466,26,6.73e-05,4.89e+3],
	    "Z" => [516,21,4.61e-05,4.10e+3],
	    "V" => [544,26,3.78e-5,3.74e+3],
	    "S" => [655,22,1.81e-05,2.60e+3],
	    );


my %SystemHash =
    (
     "Johnson" => \%Johnson,
     "johnson" => \%Johnson,

     "Cousins" => \%Cousins,
     "cousins" => \%Cousins,

     "Geneva" => \%Geneva,
     "geneva" => \%Geneva,

     "Stromgren" => \%Stromgren,
     "stromgren" => \%Stromgren,
     "uvby" => \%Stromgren,

     "2Mass" => \%TwoMass,
     "2MASS" => \%TwoMass,
     "TwoMass" => \%TwoMass,

     "CIT" => \%CIT,
     "CIT/CTIO" => \%CIT,
     "Elias" => \%CIT,

     "Walraven" => \%Walraven,
     "walraven" => \%Walraven,
     "WBVR" => \%Walraven,
     "wbvr" => \%Walraven,

     "Tycho" => \%Tycho,
     "Hipparcos" => \%Hipparcos,

     "ThuanGunn" => \%ThuanGunn,
     "Gunn" => \%ThuanGunn,

     "DDO" => \%DDO,

     "Vilnius" => \%Vilnius,

     "Straizys" => \%Straizys,

     "Eggen102" => \%Eggen102,

     "KronComet" => \%KronComet,
	 
     "Oja" => \%Oja,

     "13color" => \%thirteencolor,
	 
     "SDSS" => \%SDSS,

     );

my @SystemList = ("Johnson","Cousins","Geneva","Stromgren","2Mass",
	       "CIT","Walraven","Hipparcos","Tycho","Gunn","DDO","Vilnius","Straizys","Eggen102","KronComet","Oja","13color","SDSS");



#########################################################################
## fLambda -- compute fLambda for input photometry
#########################################################################
sub fLambda ($$$;$)
{
    my ($systemName,$band,$mag,$magError) = @_;
    $mag = numberOnly($mag);
    $mag = 0 if (!isNumeric($mag));
    my ($flux, $fluxError);
    if (!defined (my $system = $SystemHash{$systemName}))
    {
	print STDERR "Error: $systemName photometric system not defined\nSupported systems: @SystemList\n";
    }
    elsif (!defined $system->{$band})
    {
	my @tmp = keys %{ $system };
	print STDERR "Error: band $band not defined in $systemName photometric system\nSupported bands in $systemName system: @tmp\n";
    }
    else
    {
	$flux = ${ $system->{$band} }[$zeroFLambda];
	$flux *= 10**(-$mag / 2.5);
	if (isNumeric($magError))
	{
	    $fluxError = $flux * (1.0 - 10.0**(-$magError / 2.5));
	}
    }
    return ($flux,$fluxError);
} ## end fLambda  #######################################################




#########################################################################
## fNu -- compute fNu for input photometry
#########################################################################
sub fNu
{
## declarations fNu  ####################################################
    my ($systemName,$band,$mag,$magError) = @_;

## begin fNu  ###########################################################
    if (defined (my $system = $SystemHash{$systemName}))
    {
	if (grep m/$band/, keys %{ $system })
	{
	    if ((defined $mag) && ($mag =~ m/\d+/))
	    {
		my $flux = ${ $system->{$band} }[$zeroFNu]; 
		$flux *= 10**(-$mag / 2.5);
		if ((defined $magError) && ($magError =~ m/\d+/))
		{
		    my $fluxError = $flux * (1.0 - 10.0**(-$magError / 2.5));
		    return ($flux,$fluxError);
		}
		else
		{
		    return $flux;
		}
	    } ## end defined $mag
	} ## band defined
	else
	{
	    my @tmp = keys %{ $system };
	    print STDERR "Error: band $band not defined in $systemName photometric system\n";
	    print STDERR "Supported bands in $systemName system: @tmp\n";
	    return undef;
	}
    } ## end system defined
    else
    {
	print STDERR "Error: $systemName photometric system not defined\n";
	print STDERR "Supported systems: @SystemList\n";
	return undef;
    }

} ## end fNu  ###########################################################




#########################################################################
## bandParameters -- return the band parameters for a given band
#########################################################################
sub bandParameters
{
## declarations bandParameters  #########################################
    my ($systemName,$band) = @_;

## begin bandParameters  ################################################
    if (defined (my $system = $SystemHash{$systemName}))
    {
	if (grep m/$band/, keys %{ $system })
	{
	    return @{ $system->{$band} };
	} ## band defined
	else
	{
	    my @tmp = keys %{ $system };
	    print STDERR "Error: band $band not defined in $systemName photometric system\n";
	    print STDERR "Supported bands in $systemName system: @tmp\n";
	    return undef;
	}
    } ## end system defined
    else
    {
	print STDERR "Error: $systemName photometric system not defined\n";
	return undef;
    }

} ## end bandParameters  ################################################



#########################################################################
## matchBand -- try to match a photometric band from input wavelengths
#########################################################################
sub matchBand
{
## declarations matchBand  ##############################################
    my $lambda = shift @_;
    my $deltaLambda = shift @_;

    my $tol = 15.0;  ## check tolerance

## begin matchBand  #####################################################

    ## Nothing really to do except an exhausting linear search
    foreach my $system (@SystemList)
    {
	foreach my $band (keys %{ $SystemHash{$system} })
	{
	    my ($testLambda,$testDLambda,$z1,$z2)
		= &bandParameters($system,$band);

	    if ((abs($lambda - $testLambda) < $tol)
		&& (abs($deltaLambda - $testDLambda) < $tol))
	    {
		return ($system,$band);
	    } ## end returning match

	} ## end foreach $band

    } ## end foreach $system

    ## If we're still here, return undef to signal no match
    return undef;

} ## end matchBand  #####################################################





#########################################################################
## fLambdaToFNu -- convert fLambda into fNu
#########################################################################
sub fLambdaToFNu
{
## declarations fLambdaToFNu  ###########################################
    my $lambda = (shift @_) * 1.0e-9;  ## input lambda in nm, conv to m

    ## Input fLambda units are in erg / s / cm^2 / um
    ## Convert to MKS W / m^2 / m
    my $fLambda = (shift @_) * (1e-7 * 1e4 * 1e6);
    my $efLambda = shift @_;
    $efLambda *= (1e-7 * 1e4 * 1e6) if (defined $efLambda);

    my $cMKS = 2.9979e8; ## MKS speed of light (m/s)

## begin fLambdaToFNu  ##################################################

    ## Now into MKS Fnu W / m^2 / Hz
    ## fNu = fLambda * lambda^2 / c
    my $fNu = $fLambda * ($lambda**2.0) / $cMKS;

    ## Finally into Jy
    $fNu *= 1.0e26;

    if (defined $efLambda)  ## If we also have a flux error
    {
	my $efNu = $fNu * ($efLambda / $fLambda);
	return ($fNu,$efNu);  ## return both
    }

    return $fNu;  ## Otherwise just return the $fNu
	
} ## end fLambdaToFNu  ##################################################



#########################################################################
## fNuToFLambda -- convert fNu into fLambda
#########################################################################
sub fNuToFLambda
{
## declarations fNuToFLambda  ###########################################
    my $lambda = (shift @_) * 1.0e-7;  ## input lambda in nm, conv to cm

    ## Input fNu units are in Jy (10^-26 W / m^2 / Hz)
    ## Convert to CGS erg / s / cm^2 / Hz
    my $fNu = (shift @_) * (1e-26 * 1e7 / 1e4);
    my $efNu = shift @_;
    $efNu *= (1e-26 * 1e7 / 1e4) if (defined $efNu);

    my $cCGS = 2.9979e10; ## CGS speed of light (cm/s)

## begin fNuToFLambda  ##################################################

    ## Now into CGS Flambda erg / s / cm^2 / cm
    ## fLambda = fNu * c / lambda^2
    my $fLambda = $fNu * $cCGS / ($lambda**2.0);

    ## Finally into my prefered erg / s / cm^2 / um
    $fLambda *= 1e-4;

    if (defined $efNu)  ## If we also have a flux error
    {
	my $efLambda = $fLambda * ($efNu / $fNu);
	return ($fLambda,$efLambda);  ## return both
    }

    return $fLambda;  ## Otherwise just return the $fLambda
	
} ## end fNuToFLambda  ##################################################



#########################################################################
## formatPhotRec -- create a formatted photometry/flux record
## 	for use by fbol
## units 0=mag (default), 1=FLambda, 2=FNu
#########################################################################
sub formatPhotRec
{
## declarations formatPhotRec  #######################################
    my ($name,$systemName,$band,$mag,$in_err,$comment,$units) = @_;

    if (!defined $units)
    {
	$units = UNITS_MAG;
    }

    my $outRec;
## begin formatPhotRec  ##############################################
    my ($lambda,$pb,$z1,$z2) = &bandParameters($systemName,$band);
    if ($units == UNITS_FLAMBDA)  ## reporting Flambda (erg cm^-2 s^-1 um^-1)
    {
        my ($flux,$err) = &photometry::fLambda($systemName,$band,$mag,$in_err);
        my ($fluxfmt, $errfmt) = ("%.3g", "%.2g");
        ($flux, $fluxfmt) = ("?", "%s") if (!isNumeric($flux));
        ($err, $errfmt) = ("?", "%s") if (!isNumeric($err));
        $outRec = "Fl " . formatValues("?", " ", "%s", $name, "%4.0f", $lambda,
                                  "%3.0f", $pb, $fluxfmt, $flux, $errfmt, $err,
                                  "# %s", $comment) . "\n";
    }
    elsif ($units == UNITS_FNU)  ## Reporting Fnu (Jy)
    {
        my ($flux,$err) = &photometry::fNu($systemName,$band,$mag,$in_err);
        $outRec = sprintf "Fn %s %4.0f %3.0f %g %g  # $comment\n",
            $name,$lambda,$pb,$flux,$err;
    }
    else ## default -- reporting Mags
    {
        $outRec = sprintf "M %s %s %s %s %s  # %s\n", nonNull($name, "?"), nonNull($systemName, "?"), nonNull($band, "?"), nonNull($mag, "?"), nonNull($in_err, "?"), nonNull($comment, "?");
    }

    return $outRec;
} ## end formatPhotRec  ##############################################
