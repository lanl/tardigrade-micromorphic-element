class Configuration(object):
    """A class which defines the configuration used for the setup"""

    def __init__(self):
        
        ### Define the path to the eigen library ###
        self.eigen    = "/path/to/eigen"
        
        #Uncomment for Nathan Miller's laptop
        #self.eigen    = "g++"

        #Uncomment for nami2227@soils.colorado.edu location
        #self.eigen    = "/data/home/students/nami2227/Code/Eigen/eigen3.3.4"
        
        ### Define the path to the compiler library ###
        self.compiler = "/path/to/compiler"
        
        #Uncomment for Nathan Miller's laptop

        #Uncomment for nami2227@soils.colorado.edu location
        #self.compiler = "/data/home/students/nami2227/Code/gcc/gcc-7.2.0/bin/g++"
