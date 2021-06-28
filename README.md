A C++ particle physics project to implement the resonant and rare mode
Decay Model described in Chapter 6 of my thesis, https://cds.cern.ch/record/2756320?ln=en

The idea behind the project is to reweight Monte Carlo samples describing particle decays 
to include additional contributions from resonances not included within the vanilla simulation.

This is accomplished by calculating the helicity angles of the different particle daughters 
and then building an angular model to incorporate the different resonances.

This project uses the ROOT framework https://root.cern.ch/ and boost libraries https://www.boost.org/

Boost is used to read in the different config json files which store the relevant particle and
resonance information.

With thanks to Dr Tom Blake, my PhD supervisor for all of his help and work with me on this project
