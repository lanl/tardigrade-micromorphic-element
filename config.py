class Configuration(object):
    """A class which defines the configuration used for the setup"""

    def __init__(self):
        
        ### Define the path to the eigen library ###
        self.eigen    = "/path/to/eigen"
        
        #Uncomment for Nathan Miller's laptop
        #self.eigen    = "C:\usr\share\cpp\eigen3.3.4"

        #Uncomment for nami2227@soils.colorado.edu location
        #self.eigen    = "/data/home/students/nami2227/Code/Eigen/eigen3.3.4"
        
        ### Define the path to the compiler library ###
        self.compiler = "/path/to/compiler"
        
        #Uncomment for Nathan Miller's laptop
        #self.compiler = "g++"

        #Uncomment for nami2227@soils.colorado.edu location
        #self.compiler = "/data/home/students/nami2227/Code/gcc/gcc-7.2.0/bin/g++"

        ### Define the make command
        self.make = "make"
        
        #Uncomment for Nathan Miller's laptop
        #self.make = "mingw32-make"

    def __repr__(self):
        output_string    = "\n######################################\n"+\
                           "###     Configuration Settings     ###\n"+\
                           "######################################\n\n"+\
                           "Configuration settings for the\n"+\
                           "required settings for the micromorphic\n"+\
                           "finite element.\n"
        eigen_string     = "Eigen Library: " + str(self.eigen) + "\n"
        compiler_string  = "Compiler Command: " + str(self.compiler) + "\n"
        make_string      = "Make Command: " + str(self.make) + "\n"

        output_string   += "\n### Libraries ###\n\n"+eigen_string+\
                           "\n### Compiler Settings ###\n\n"+compiler_string+\
                           "\n### Make Settings ###\n\n"+make_string
        return output_string

if __name__=="__main__":
    print Configuration()
