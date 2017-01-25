# muon_inflight_decay
A simulation that will simulate cosmic muon in flight decays as they travel through ice.

This project was my senior thesis at Université Libre de Bruxelles working under the lead data acquisition scientist 
for the ICECUBE collaboration. ICECUBE is a collection of equally spaced detectors occupying a one cubic kilometer volume
in the antarctic ice. Neutrinos can't be detected directly, but when a neutrino passes through a dielectric medium, it can
interact with the atomic nuclei and produce a charged lepton which emits Cherenkov radiation if the neutrino is sufficiently
high energy. However, muons arising as the product of interaction from high energy cosmic particles with the atmosphere reach
the ground with speeds quite close to light and decay into an electron, an electron antineutrino, and a muon neutrino
(antimuons decay into the corresponding antiparticles). These decay events will also produce Cherenkov radiation
and so they must be accounted for to distinguish neutrino candidates.

Without going into too much detail, the lifetime of the muon from the reference frame of the laboratory, in this case
the detector array, is dependent on its speed, which in turn depends on its energy. As the muon travels through ice,
it loses energy by interacting with the water molecules, therefore reducing its speed and subsequently changing its
probability to decay. 

This simulation is designed to model the distribution of muons and their energies that would decay inside ICECUBE.
