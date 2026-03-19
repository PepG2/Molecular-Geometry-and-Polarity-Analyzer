Molecular Geometry and Polarity Analyzer 

By Ethan Chen



  Molecular geometry and polarity determine many properties in molecules such as melting points, boiling points, crystal shape, solubility, intermolecular forces, and alignment in electrical fields. Traditional polarity predictions are based on qualitative symmetry reasoning and may not work well for more complex molecules. This project aims to tackle this problem through a quantitative stance, it predicts and draws molecular geometries with the use of vectors allowing for calculations of molecular polarities through geometry and bond properties.

  VSEPR theory shows how electrons arrange themselves in bonds and lone pairs in order to stay as far away from other electrons as possible. Using this we were able to define a 
repulsion energy and use force to rotate the electron domains iteratively to minimize repulsion and create ideal geometries. Lone pairs repel more than bonds so we used weights to model this. Furthermore, bonds have electronegativity differences that determine the bond's polarity and so using this and the ideal geometry, we are able to determine bond dipole moments and therefore; the molecular dipole moment. 

  Starting with a number of random vectors representing bonds determined from bonding theory, we are able to define a repulsion energy for one bond by summing the reciprocals of the distances between it and other electron domains that are multiplied by lone pair and bond interaction weights obtained through calibration using known experimental bond angles. Using force, the negative gradient of energy, we are able to model in which direction the repulsion would push on our bond. By making thousands of little rotations that maximize repulsion, we are eventually able to create an ideal geometry when the repulsion energy is essentially no longer changing and the forces are essentially 0. Using the bonds generated from our created geometry and the electronegativity difference of each individual bond, we can find each bond dipole moment. Finally, we can then proceed to add the dipole moments to find the molecular dipole moment. 

