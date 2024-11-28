========================================================================
    STATIC LIBRARY : atCLAMLink Project Overview
========================================================================

AppWizard has created this atCLAMLink library project for you. 
No source files were created as part of your project.


atCLAMLink.vcproj
    This is the main project file for VC++ projects generated using an Application Wizard. 
    It contains information about the version of Visual C++ that generated the file, and 
    information about the platforms, configurations, and project features selected with the
    Application Wizard.

/////////////////////////////////////////////////////////////////////////////
Other notes:

AppWizard uses "TODO:" comments to indicate parts of the source code you
should add to or customize.

/////////////////////////////////////////////////////////////////////////////


NOTES:
The requirements of this module are at https://wiki.csiro.au/confluence/display/ML/Proof+of+Concept+Requirements.



ok.
it would be good ot have a mathematical way of getting the 'best' outcome.
If the indicators is greater than the ideal value then we want to select a strategy that will have the effect of decreasing the indicator.
There should be a decision about if we should be trying to get a large decrease or a small decrease as well. This will depend on the how much the
indicator is outside the 'good' range.

We could calculate the percentage that the indicator is out of the range.

For example for salinity the acceptable indicator range might be 
34-36

If the salinity got to 36.5 then this corresponds to a desired decrease of:

desired decrease 
= ((36-34)/2 - 36.5)/36.5 
= (35-36.5)/36.5
= 2.5/36.5
= -0.07%

So a small decrease of 10-20 percent would probably be perfect 
but a large increase would push it back under the acceptable range and we might end up with some instability.




so we want to find the strategy that has the highest probability to creating this desired change.

so we need to assign a weight to each of the impacts that is dependant on the desired change.

So.
Calculate the desired change:
Define impact value 

So if we are after a change of -7% then the impact levels would be:

Large Decrease 20-50% 1
Small Decrease 10-20% 2
No Change 0
Small increase 10-20% -1
Large increase 20-50% -2

We might need to play around with the weightings. Some managers might be conservative and assigned a higher weight to no change 
instead of pushing for a larger change than is required.

Then select the strategy where:

the summation is the highest.

So biggest challenge is working out the different weightings to assign to the possible impact values.


Not sure what to do about the objective weighting. These don't really seem to have much relevance to the indicator driven mangement.


Another question is if we want the selection of strategy to be random or to set up an order of things that we try.
We almost want a way of setting up a weighting of indicator, mangement scenario and strategy.
Or at last indicator and management scenario. Otherwise we might be trying strategies than could not possibly have a effect on the indicator.

 
We will have to store additional information about the output states in the code. The possible outcomes vary based on which output option is selected.
Refer to page 26 of the CLAM report for more information.
This will make it harder to select the correct weight for the output.

OK.
Indicator information is read in .
Now need to link the indicators to the clam scenario. For now i might just link the indicators with a single sceanrio and run all options in that scenario and select the one that 
has the 'best' outcome.



TODO;

1. Link indicators with output values from the CLAM - should perhaps read this in from the input file. - DONE
- Need to setup ranges for each output state so we can work out the weighting of the impacts based on the current required change.
Should deal with the simple case where there are set percentage changes for now. Can work on the standard dev output values later as this is more complicated.

2. Need to write code to assign the impact values based on required change and output states 
2. Write code to generate the CLAM input file based on select scenarios
