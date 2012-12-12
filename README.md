SEDfitspTypeFITS
================

Credit: Description of the files in this repository is credited to Gerard van Belle of the Lowell Observatory.

13-color_source_120806a.txt and MSX-source-120810a.txt - Codes that collects the photometry from the GCPD that the 
file photometry.batx could not (The GCPD lookup could not work for the former; and the VizieR lookup did not work 
for the latter.)

GCPD2.py - A python script that queries the GCPD database for known photometric surveys and their results.

getAllPhotometry.batx - This is a simple shell script that calls the previously installed fbol package to query 
the SIMBAD database for photometry, also calls GCPD2.py.

photometry.pm - Perl script (application/x-perl) that defines useful photometric quantities such as wavelengths, 
passbands, zeropoints and the like of.
