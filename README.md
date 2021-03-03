# eye_processing
Processing of eye tracking data to obtain various behavioral measures of interest (looking time, saccadic latencies, habituation time).

Integrating data collected from eye tracking with data coded on boris (https://www.boris.unito.it/)

Analyzing what participants have seen with an information-theory based computational model (ITDmodel.m)

The input is eye-tracking data pre-processed using the fixation-extraction algorithm I2MC (https://github.com/royhessels/I2MC)

- Adult_processing.m is the main script to analyze data from adults
- Infant_processing.m is the main script to analyze data from infants
- boris2R.m and postboris.m integrate data coded on boris to the infant script

