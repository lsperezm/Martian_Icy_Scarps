Started with v12 from the files uploaded to Figshare for Bramson et al. 2019

-----

Started with v10_z_A6-atmos-addsbackinRerad but made major modifications so that it works for thin lags

Gets rid of extraneous, old parts of code too.

This version combines v11 with Free/Forced Convection basic thermal model to calculate the flat, regional values first
and then calculates and saves surface sublimation at each orbital timestep. Makes sure to reorder arrays to be ordered
in diurnal loops for the surface sublimation calculations.

This version calls several of the scripts to now be in functions, which are supposed to run faster, and will also
make sure nothing is inadvertantly carried over between timesteps.

It also fixes some issues in calcOrbitalParams where there was the possibility of getting imaginary numbers
in the solar azimuth angle which could then propogate into temperatures and all other calculations. It adds
checks to cosi and cos of the azimuth angle to make sure they stay between [-1,1].

As setup now with output, V12 is best used for just calculating free and forced sublimation values throughout
orbital histories. V13 will have a more comprehensive output that compares sublimation through a diffusive lag
to the outputs from runs of this V12 code (or maybe I'll decide to have it calculate the surface sublimation
values so the code is more comprehensive, even though it will make it slower).

MAY 29, 2018:
UPDATED TO MULTIPLY BY EarthYearLength/MarsYearLength so that the annual result is in per Earth Years not Mars Years,
as to use with Laskar timesteps and compare to the diffusive ice loss