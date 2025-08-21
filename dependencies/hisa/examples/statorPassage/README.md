# Body Force Model example case

This is the verification case for a Body Force Model implementation.

## Body Force Model (BFM)
The current implementation is the Hall-Thollet BFM, which is the 
extension of Hall's BFM to compressible flows and added metal blockage terms 
effects by Thollet.

Sources:

Hall, David Kenneth. *Analysis of civil aircraft propulsors with boundary layer 
ingestion*. Diss. Massachusetts Institute of Technology, 2015. [Link](https://dspace.mit.edu/bitstream/handle/1721.1/97353/910627725-MIT.pdf?sequence=1&isAllowed=y)

Thollet, William. *Body force modeling of fan-airframe interactions*. Diss. 
ISAE-SUPAERO, 2017. [Link](https://depozit.isae.fr/theses/2017/2017_Thollet_William_A.pdf)


## Example case
The original BFM verification case (including the sym_stator_BFM.vtk mesh file) was kindly provided by 
Evert Bunschoten, who employed the same verification case for verifying a BFM 
implementation in SU2. 

[Gitlab repository](https://github.com/EvertBunschoten/Meangen2BFM)

[Documentation](https://repository.tudelft.nl/islandora/object/uuid:eb03a5d5-915c-47f8-a586-f9d5b175c6ba/datastream/OBJ/download)
