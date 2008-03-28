  pbflags.MonitorCorrect=0
This flag is to perform corrections for Monitor before the front 3He Cell.
The default is 0
   pbflags.PolMonitorCorrect=1
This flag is to perform corrections for Monitor after the front 3He Cell.
The default is 1.  If it is 0 and the MonitorCorrect flag is also 0, then the program will assume that you are
counting by time
   pbflags.MonoSelect=1
Assumes that there is a PG filter before the front 3He Cell and uses that to
calculate the amount of lambda/2 in the beam.  If it is set to 0, there is no
filter and a correction is performed based on an empirical fit for the amount of
lambda/2 in the beam
   pbflags.NoNegativeCS=0
Use this flag to constrain all cross sections to be non negative.  Hopefully not
necessary
For the following the order of the flags is:
pp mm pm mp
   pbflags.CountsEnable=[0,0,0,0]
A one will enable the given channel for correction. The program will expect a
file to be read for this data
   pbflags.CountsAdd1=[0,0,0,0]
Only the first 2 indices of this flag are used
This is used to add cross sections together to get better statistics
For example:
CountsAdd1 1 2 0 0
will take the data from the -- cross section and add it to the ++ cross section
and use that sum for corrections involving the ++ channel (that is data in the
first channel is replaced)
   pbflags.CountsAdd2=[0,0,0,0]
?
   pbflags.Sconstrain=[0,0,0,0]
Use this to constrain a channel
   pbflags.Spp=[0,0,0,0]
Constraint coefficients for the ++ channel
For example:
Sconstrain 1 0 0 0
Spp 0 1 0 0
Will constrain S++=S--
Sconstrain 0 1 0 0
Smm 1 0 0 0
Will constrain S--=S++

Sconstrain 1 1 0 0
without any flags set will simply set S++=S--=0
etc.

   pbflags.Smm=[0,0,0,0]
Constraint coefficients for the -- channel
   pbflags.Spm=[0,0,0,0]
Constraint coefficients for the +- channel
   pbflags.Smp=[0,0,0,0]
Constraint coefficients for the -+ channel




----- End forwarded message -----


