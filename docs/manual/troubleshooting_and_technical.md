Troubleshooting
===============

Logging
-------
More verbose logging can be provided with an EasyLogging log config file.
The file may look like this:

	* GLOBAL:
	   FORMAT               =  "%datetime %msg"
	   FILENAME             =  "cagee.log"
	   ENABLED              =  true
	   TO_FILE              =  true
	   TO_STANDARD_OUTPUT   =  true
	   SUBSECOND_PRECISION  =  6
	   PERFORMANCE_TRACKING =  true
	   MAX_LOG_FILE_SIZE    =  2097152 ## 2MB - Comment starts with two hashes (##)
	   LOG_FLUSH_THRESHOLD  =  100 ## Flush after every 100 logs
	* DEBUG:
	   FORMAT               = "%datetime{%d/%M} %func %msg"

For more information, see https://github.com/amrayn/easyloggingpp#using-configuration-file

Pass the config file to CAGEE with the --log_config flag. For example,

    cagee -c cagee_estimate.cfg --log_config log.config	   

Technical
=========

How does the optimizer work?
----------------------------

The Nelder-Mead optimization algorithm is used. It runs until it can
find a difference of less than 1e-6 in either the calculated score or
the calculated value, or for 10,000 iterations. The parameters that are
used for the optimizer are as follows:

-   rho: 1 (reflection)

-   chi: 2 (expansion)

-   psi: 0.5 (contraction)

-   sigma: 0.5 (shrink)

In some cases, the optimizer suggests values that cannot be calculated
(due to saturation, negative values, or other reasons) In this case, an
infinite score is returned and the optimizer continues.
