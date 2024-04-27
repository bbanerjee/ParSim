#/bin/env python
import numpy as np
import xml.etree.ElementTree as ET
import argparse
import os,sys
import re,operator
import time
import cStringIO

try:
    # Import mpi4py; if it fails, execute serially without it
    from mpi4py import MPI
except ImportError:
    hasmpi = False
    print "WARNING: mpi4py could not be imported.  Using single rank execution."
else:
    hasmpi = True

"""
Filename:   uda2vtk.py

Summary:    Conversion utility to read Uintah *.uda simulation output files and 
            create VTK files for grid, particle and 
            particle domain visualization and analysis

Usage:      python uda2vtk.py </path/filename.uda> [-flags]

Output:     VTK XML files created for reading with ParaView or other visualization tools 
            grid.pvd                    Uintah traditional grid files 
            particles.pvd               Uintah traditional particle centroid files 
            domains.pvd                 Uintah particle domain files 
            
Revisions:  YYMMDD    Author            Comments
-------------------------------------------------------------------------------------
            140505    Brian Leavy       created uda2vtk.py
            130601    John Schmidt      started Uintah uda reader in python 
            141007    Sean Ziegler      mpi4py processing of timesteps 
"""
# Constants
tolerance = 1.e-100                                                             # set numbers smaller than this to zero
scalar = 1                                                                      # length of lists for different types   
vector = 3
tensor = 9

# Utility functions
def chunks(lst,n):
    """ Partitions string list into n-sized chunks of floats """
    for i in xrange(0,len(lst),n):
      yield map(float,lst[i:i+n])

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def IntVector(iv):
    s = iv[1:-1]
    s_array = s.split(',')
    int_vector = map(int,s_array)
    return int_vector

def Vector(iv):
    s = iv[1:-1]
    s_array = s.split(',')
    vec = map(float,s_array)
    return vec

def find_type(type):
    m = re.search("<[a-zA-Z]+[\d]*>",type)
    t = m.group()
    return t[1:-1]

def generate_coordinates(high_index,low_index,upper,lower):
    num_x = high_index[0] - low_index[0]
    num_y = high_index[1] - low_index[1]
    num_z = high_index[2] - low_index[2]
    space_x = (upper[0] - lower[0])/num_x
    space_y = (upper[1] - lower[1])/num_y
    space_z = (upper[2] - lower[2])/num_z
    x_coord = []
    y_coord = []
    z_coord = []
    for i in range(num_x+1):
        x_coord.append(lower[0] + float(i)*space_x)
    for i in range(num_y+1):
        y_coord.append(lower[1] + float(i)*space_y)
    for i in range(num_z+1):
        z_coord.append(lower[2] + float(i)*space_z)

    return (x_coord,y_coord,z_coord)


def create_vtkfile_element(vtk_type,version=None,endianness=None):
    vtkfile_element = ET.Element('VTKFile')
    vtkfile_element.attrib['type'] = vtk_type
    if version is not None:
        vtkfile_element.attrib['version'] = version
    else:
        # VTK >=4: vtkfile_element.attrib['version'] = '0.1'
        # VTK >=6: vtkfile_element.attrib['version'] = '1.0'
        vtkfile_element.attrib['version'] = '1.0'

    if endianness is not None:
        vtkfile_element.attrib['byte_order'] = endianness
    else:
        vtkfile_element.attrib['byte_order'] = 'LittleEndian'
        
    return vtkfile_element



def calcDomain(xp,r1,r2,r3,interpolator):
    """
    Calculate particle domain corners 
    INPUT:    xp                        particle centroid 
              r1, r2, r3                current Rvectors 
              interpolator              interpolator type 
    OUTPUT:   xpc                       particle corners 
    """ 
    if interpolator == 'cpti':                                                  # CPTI
      cpp=4                                                                     # corners per particle 
    else:                                                                       # CPDI, CBDI, GIMP
      cpp=8

    dim=len(xp)                                                                 # dimensions
    xpc=np.zeros((cpp,dim),float)                                               # particle domain corners
    # MODIFIED VERSION ALREADY HAS APPLIED F, Rvectors, spacing?
    # deformation gradient
    # Convected Particle Domain (CPDI or Tetrahedral/Triangular CPTI) Interpolation
    if interpolator == 'cpti':                                                  # CPTI
      # calculate the coordinates of 4 corners of the particle domain       
      xpc[0]=(xp[0]-(r1[0]+r2[0]+r3[0])/4., xp[1]-(r1[1]+r2[1]+r3[1])/4., xp[2]-(r1[2]+r2[2]+r3[2])/4.)
      xpc[1]=(xpc[0,0]+r1[0], xpc[0,1]+r1[1], xpc[0,2]+r1[2])
      xpc[2]=(xpc[0,0]+r2[0], xpc[0,1]+r2[1], xpc[0,2]+r2[2])
      xpc[3]=(xpc[0,0]+r3[0], xpc[0,1]+r3[1], xpc[0,2]+r3[2])
    else:                                                                       # CPDI, CBDI, GIMP
      # calculate the coordinates of 8 corners of the particle domain
      xpc[0]=(xp[0]-r1[0]/2.-r2[0]/2.-r3[0]/2., xp[1]-r1[1]/2.-r2[1]/2.-r3[1]/2., xp[2]-r1[2]/2.-r2[2]/2.-r3[2]/2.)
      xpc[1]=(xp[0]+r1[0]/2.-r2[0]/2.-r3[0]/2., xp[1]+r1[1]/2.-r2[1]/2.-r3[1]/2., xp[2]+r1[2]/2.-r2[2]/2.-r3[2]/2.)
      xpc[2]=(xp[0]+r1[0]/2.+r2[0]/2.-r3[0]/2., xp[1]+r1[1]/2.+r2[1]/2.-r3[1]/2., xp[2]+r1[2]/2.+r2[2]/2.-r3[2]/2.)
      xpc[3]=(xp[0]-r1[0]/2.+r2[0]/2.-r3[0]/2., xp[1]-r1[1]/2.+r2[1]/2.-r3[1]/2., xp[2]-r1[2]/2.+r2[2]/2.-r3[2]/2.)
      xpc[4]=(xp[0]-r1[0]/2.-r2[0]/2.+r3[0]/2., xp[1]-r1[1]/2.-r2[1]/2.+r3[1]/2., xp[2]-r1[2]/2.-r2[2]/2.+r3[2]/2.)
      xpc[5]=(xp[0]+r1[0]/2.-r2[0]/2.+r3[0]/2., xp[1]+r1[1]/2.-r2[1]/2.+r3[1]/2., xp[2]+r1[2]/2.-r2[2]/2.+r3[2]/2.)
      xpc[6]=(xp[0]+r1[0]/2.+r2[0]/2.+r3[0]/2., xp[1]+r1[1]/2.+r2[1]/2.+r3[1]/2., xp[2]+r1[2]/2.+r2[2]/2.+r3[2]/2.)
      xpc[7]=(xp[0]-r1[0]/2.+r2[0]/2.+r3[0]/2., xp[1]-r1[1]/2.+r2[1]/2.+r3[1]/2., xp[2]-r1[2]/2.+r2[2]/2.+r3[2]/2.)
      # calculate the coordinates of 8 corners of the particle domain from lower corner
      #xpc[0]=(xp[0], xp[1], xp[2])
      #xpc[1]=(xp[0]+r1[0], xp[1]+r1[1], xp[2]+r1[2])
      #xpc[2]=(xp[0]+r1[0]+r2[0], xp[1]+r1[1]+r2[1], xp[2]+r1[2]+r2[2])
      #xpc[3]=(xp[0]+r2[0], xp[1]+r2[1], xp[2]+r2[2])
      #xpc[4]=(xp[0]+r3[0], xp[1]+r3[1], xp[2]+r3[2])
      #xpc[5]=(xp[0]+r1[0]+r3[0], xp[1]+r1[1]+r3[1], xp[2]+r1[2]+r3[2])
      #xpc[6]=(xp[0]+r1[0]+r2[0]+r3[0], xp[1]+r1[1]+r2[1]+r3[1], xp[2]+r1[2]+r2[2]+r3[2])
      #xpc[7]=(xp[0]+r2[0]+r3[0], xp[1]+r2[1]+r3[1], xp[2]+r2[2]+r3[2])

    return xpc

def oldCalcDomain(xp,defgrad,oldr1,oldr2,oldr3,interpolator):
    """
    Calculate particle domain corners 
    INPUT:    xp                        particle centroid 
              oldr1, oldr2, oldr3       current Rvectors 
              defgrad, interpolator     deformation gradient and interpolator type 
    OUTPUT:   xpc                       particle corners 
    """ 
    if interpolator == 'cpti':                                                  # CPTI
      cpp=4                                                                     # corners per particle 
    else:                                                                       # CPDI, CBDI, GIMP
      cpp=8

    dim=len(xp)                                                                 # dimensions
    xpc=np.zeros((cpp,dim),float)                                               # particle domain corners
    # MODIFIED VERSION ALREADY HAS APPLIED F, Rvectors, spacing?
    # deformation gradient
    F=np.eye(dim)
    F[0]=defgrad[0],defgrad[1],defgrad[2] 
    F[1]=defgrad[3],defgrad[4],defgrad[5] 
    F[2]=defgrad[6],defgrad[7],defgrad[8] 
    # updated rvectors based on deformation gradient
    r1=np.dot(F,oldr1) 
    r2=np.dot(F,oldr2) 
    r3=np.dot(F,oldr3)
    # Convected Particle Domain (CPDI or Tetrahedral/Triangular CPTI) Interpolation
    if interpolator == 'cpti':                                                  # CPTI
      # calculate the coordinates of 4 corners of the particle domain       
      xpc[0]=(xp[0]-(r1[0]+r2[0]+r3[0])/4., xp[1]-(r1[1]+r2[1]+r3[1])/4., xp[2]-(r1[2]+r2[2]+r3[2])/4.)
      xpc[1]=(xpc[0,0]+r1[0], xpc[0,1]+r1[1], xpc[0,2]+r1[2])
      xpc[2]=(xpc[0,0]+r2[0], xpc[0,1]+r2[1], xpc[0,2]+r2[2])
      xpc[3]=(xpc[0,0]+r3[0], xpc[0,1]+r3[1], xpc[0,2]+r3[2])
    else:                                                                       # CPDI, CBDI, GIMP
      # calculate the coordinates of 8 corners of the particle domain
      xpc[0]=(xp[0]-r1[0]/2.-r2[0]/2.-r3[0]/2., xp[1]-r1[1]/2.-r2[1]/2.-r3[1]/2., xp[2]-r1[2]/2.-r2[2]/2.-r3[2]/2.)
      xpc[1]=(xp[0]+r1[0]/2.-r2[0]/2.-r3[0]/2., xp[1]+r1[1]/2.-r2[1]/2.-r3[1]/2., xp[2]+r1[2]/2.-r2[2]/2.-r3[2]/2.)
      xpc[2]=(xp[0]+r1[0]/2.+r2[0]/2.-r3[0]/2., xp[1]+r1[1]/2.+r2[1]/2.-r3[1]/2., xp[2]+r1[2]/2.+r2[2]/2.-r3[2]/2.)
      xpc[3]=(xp[0]-r1[0]/2.+r2[0]/2.-r3[0]/2., xp[1]-r1[1]/2.+r2[1]/2.-r3[1]/2., xp[2]-r1[2]/2.+r2[2]/2.-r3[2]/2.)
      xpc[4]=(xp[0]-r1[0]/2.-r2[0]/2.+r3[0]/2., xp[1]-r1[1]/2.-r2[1]/2.+r3[1]/2., xp[2]-r1[2]/2.-r2[2]/2.+r3[2]/2.)
      xpc[5]=(xp[0]+r1[0]/2.-r2[0]/2.+r3[0]/2., xp[1]+r1[1]/2.-r2[1]/2.+r3[1]/2., xp[2]+r1[2]/2.-r2[2]/2.+r3[2]/2.)
      xpc[6]=(xp[0]+r1[0]/2.+r2[0]/2.+r3[0]/2., xp[1]+r1[1]/2.+r2[1]/2.+r3[1]/2., xp[2]+r1[2]/2.+r2[2]/2.+r3[2]/2.)
      xpc[7]=(xp[0]-r1[0]/2.+r2[0]/2.+r3[0]/2., xp[1]-r1[1]/2.+r2[1]/2.+r3[1]/2., xp[2]-r1[2]/2.+r2[2]/2.+r3[2]/2.)
      # calculate the coordinates of 8 corners of the particle domain from lower corner
      #xpc[0]=(xp[0], xp[1], xp[2])
      #xpc[1]=(xp[0]+r1[0], xp[1]+r1[1], xp[2]+r1[2])
      #xpc[2]=(xp[0]+r1[0]+r2[0], xp[1]+r1[1]+r2[1], xp[2]+r1[2]+r2[2])
      #xpc[3]=(xp[0]+r2[0], xp[1]+r2[1], xp[2]+r2[2])
      #xpc[4]=(xp[0]+r3[0], xp[1]+r3[1], xp[2]+r3[2])
      #xpc[5]=(xp[0]+r1[0]+r3[0], xp[1]+r1[1]+r3[1], xp[2]+r1[2]+r3[2])
      #xpc[6]=(xp[0]+r1[0]+r2[0]+r3[0], xp[1]+r1[1]+r2[1]+r3[1], xp[2]+r1[2]+r2[2]+r3[2])
      #xpc[7]=(xp[0]+r2[0]+r3[0], xp[1]+r2[1]+r3[1], xp[2]+r2[2]+r3[2])

    return xpc



# *** Uda Class ***
class Uda:
    '''Uintah uda class'''

    def __init__(self,uda=None,timedir=None,mpicomm=None):
        self.name = uda
        self.timedir = timedir
        self.mpicomm = mpicomm
        if mpicomm is not None:
          self.mpirank = mpicomm.Get_rank()
          self.mpisize = mpicomm.Get_size()
        else:
          self.mpirank, self.mpisize = 0, 1
        self.output_grid_file_names = []
        self.output_grid_mat_file_names = []
        self.output_particle_file_names = []
        self.output_particle_mat_file_names = []
        self.output_particle_domain_file_names = []
        self.output_particle_domain_mat_file_names = []
        self.materials = []
        self.gmaterials = []
        self.interpolator = ''

    def read(self,uda=None,timedir=None):
        if uda is not None:
            self.name = uda
        elif self.name is None:
            print 'Please specify an Uda directory name'
            return
            
        os.chdir(self.name)
        print '\nReading Uintah file.'
        print '  ',self.name
        try:
            self.index_tree = ET.parse("index.xml")
        except Exception:
            print 'Unexpected error opening %s' % 'index.xml' 
            return

        #print "Global variables:"
        self.read_globals()
        # Read all timesteps or specific one
        print "\nReading timesteps."
        if timedir is not None:
            self.timedir = timedir
        self.read_timesteps()

        return None

            
    def read_globals(self):
        self.global_variables = {}
        elems = self.index_tree.getroot()
        global_vars = elems.find('globals')
        global_variables_iter = global_vars.getiterator('variable')
        for global_variables in global_variables_iter:
          variable = global_variables.get('href')
          name = global_variables.get('name')
          self.global_variables[name] = variable
        return None


    def read_timesteps(self):
        self.timesteps = []
        self.timestepdirs = []
        elems = self.index_tree.getroot()
        timestep_iter = elems.getiterator('timestep') 
        i = -1
        for ts in timestep_iter:
          i += 1
          timestep_href = ts.get('href')
          timestep = timestep_href.split('/')

          # Choose a timestep for this rank based on modulus;
          #  otherwise, assign None to the timestep array element
          if i % self.mpisize == self.mpirank: 
            tstepdat = Timestep(timestep_href)
          else:
            tstepdat = None

          if self.timedir:
            if (timestep[0]==self.timedir):
              self.timestepdirs.append(timestep[0])
              self.timesteps.append(tstepdat)
              return
            else:
              pass
          else:
            self.timestepdirs.append(timestep[0])
            self.timesteps.append(tstepdat)


    def output_vtk_grid(self,filename = None):
        root_element = self.create_root_element(vtk_type = 'Collection')
        self.timestep_vtk_grid()
        self.output_vtk_grid_file(root_element)

    def output_vtk_particle(self,filename = None):
        root_element = self.create_root_element(vtk_type = 'Collection')
        self.timestep_vtk_particle()
        self.output_vtk_particle_file(root_element)

    def output_vtk_particle_domain(self,filename = None,nprocs=1):
        root_element = self.create_root_element(vtk_type = 'Collection')
        self.timestep_vtk_particle_domain(nprocs)
        self.output_vtk_particle_domain_file(root_element)

    def create_root_element(self,vtk_type):
        elems = self.index_tree.getroot()
        meta_tag = elems.find('Meta')
        endianness = meta_tag.findtext('endianness')
        if endianness == 'little_endian':
            endian_tag = 'LittleEndian'
        elif endianness == 'big_endian' :
            endian_tag = 'BigEndian'
        #root_element = create_vtkfile_element(vtk_type,'0.1',endian_tag)
        root_element = create_vtkfile_element(vtk_type,'1.0',endian_tag)

        return root_element


    # Grid vtk timesteps
    def timestep_vtk_grid(self):
        t=-1
        for ts in self.timesteps:
            t += 1
            if t % self.mpisize != self.mpirank: continue   # Only process this rank's tsteps
            tdir=self.timestepdirs[t]
            print '   grid timestep: ',tdir
            filename_grid_relative = tdir+'/grid.pvti'
            filename_grid = self.name + '/'+filename_grid_relative
            gmat_subfiles=[]
            for level in ts.grid.levels:
                os.chdir(self.name+'/'+tdir)
                patches = level.get_patches()
                # rbl: can save space if just save last material number (all materials combined):  
                #for gmat_index in ts.gmaterials:
                for patch in patches:
                  gmat_index = ts.gmaterials[-1]
                  filename_grid_mat_relative = 'l'+str(level.id)+\
                    '/grid_m'+str(gmat_index)+'_p'+str(patch.id)+'.vti'
                  filename_grid_mat = self.name+ '/' + tdir + '/'+filename_grid_mat_relative
                  self.output_grid_mat_file_names.append(filename_grid_mat_relative)
                  gmat_subfiles.append(filename_grid_mat_relative)
                  ts.output_vtk_grid(filename_grid_mat,patch,gmat_index)
            # output parallel grid.pvti for patches and all materials together 
            self.output_grid_file_names.append(filename_grid_relative)
            ts.output_vtk_pvti_grid(filename_grid,gmat_subfiles)

    # Particle centroid vtk timesteps
    def timestep_vtk_particle(self):
        t=-1
        for ts in self.timesteps:
            t += 1
            if t % self.mpisize != self.mpirank: continue   # Only process this rank's tsteps
            tdir=self.timestepdirs[t]
            print '   particle centroid timestep: ',tdir
            #filename_particle = self.name+ '/' + tdir + '/particles.pvtu'
            filename_particle_relative = tdir + '/particles.pvtu'
            filename_particle = self.name + '/'+filename_particle_relative
            mat_subfiles=[]
            for level in ts.grid.levels:
              os.chdir(self.name+'/'+tdir)
              patches = level.get_patches()
              for mat_index in ts.materials:
                for patch in patches:
                  filename_particle_mat_relative = 'l'+str(level.id)+\
                    '/particles_m'+str(mat_index)+'_p'+str(patch.id)+'.vtu'
                  filename_particle_mat = self.name+ '/' + tdir + '/'+filename_particle_mat_relative
                  self.output_particle_mat_file_names.append(filename_particle_mat_relative)
                  mat_subfiles.append(filename_particle_mat_relative)
                  ts.output_vtk_particle(filename_particle_mat,patch,mat_index)
            # output parallel particles.pvtu for all patches and materials  
            self.output_particle_file_names.append(filename_particle_relative)
            ts.output_vtk_pvtu_particle(filename_particle,mat_subfiles)


    # Particle domain vtk timesteps
    def timestep_vtk_particle_domain(self,nprocs):
        t=-1
        for ts in self.timesteps:
            t += 1
            if t % self.mpisize != self.mpirank: continue   # Only process this rank's tsteps
            #numTimesteps=len(self.timesteps)
            tdir=self.timestepdirs[t]
            print '   particle domain timestep: ',tdir
            mat_subfiles=[]
            filename_particle_domain_relative = tdir + '/domains.pvtu'
            filename_particle_domain = self.name + '/'+filename_particle_domain_relative
            for level in ts.grid.levels:
              os.chdir(self.name+'/'+tdir)
              patches = level.get_patches()
              for mat_index in ts.materials:
                for patch in patches:
                  filename_particle_domain_mat_relative = 'l'+str(level.id)+\
                    '/domains_m'+str(mat_index)+'_p'+str(patch.id)+'.vtu'
                  filename_particle_domain_mat = self.name+ '/' + tdir + '/'+filename_particle_domain_mat_relative
                  self.output_particle_domain_mat_file_names.append(filename_particle_domain_mat_relative)
                  mat_subfiles.append(filename_particle_domain_mat_relative)
                  ts.output_vtk_particle_domain(filename_particle_domain_mat,patch,mat_index)
            # Output parallel domains.pvtu for all materials
            self.output_particle_domain_file_names.append(filename_particle_domain_relative)
            ts.output_vtk_pvtu_domain(filename_particle_domain,mat_subfiles)
            
    def output_vtk_grid_file(self,root_element):
        # rbl: changes for relative paths?
        os.chdir(self.name)
        collection = ET.SubElement(root_element,"Collection")
        for i in range(len(self.output_grid_file_names)):
            data_set = ET.SubElement(collection,"DataSet")
            data_set.attrib['timestep'] = str(self.timesteps[i].currentTime)
            data_set.attrib['group'] = 'Grid elements'
            data_set.attrib['part'] = '0'
            data_set.attrib['file'] = str(self.output_grid_file_names[i])
        indent(root_element)         
        #ET.dump(root_element)
        tree = ET.ElementTree(root_element)
        # If do_compression, add? compressor="vtkZLibDataCompressor"
        #tree.write(self.name + '/grid.pvd',encoding="UTF-8",xml_declaration=True)
        tree.write('grid.pvd',encoding="UTF-8",xml_declaration=True)


    def output_vtk_particle_file(self,root_element):
        # All time steps put their valid time steps into a list of times and files
        all_times_files = []
        i = 0
        for t in range(len(self.timesteps)):
           if self.timesteps[t] is not None:
             tsd = { 't': t, 
                     'timestep': self.timesteps[t].currentTime,
                     'file': self.output_particle_file_names[i] }
             i += 1
             all_times_files.append(tsd)

        # Gather all valid timesteps from all tasks
        if self.mpicomm is not None:
          all_times_files = self.mpicomm.gather(all_times_files, root=0)

        if self.mpirank == 0: 
          # Sort all timesteps by time index into new list
          sort_times_files = [None]*len(self.timesteps)
          for rnk in all_times_files:
            for tsd in rnk:
              sort_times_files[tsd['t']] = tsd

        collection = ET.SubElement(root_element,"Collection")
        for i in range(len(self.timesteps)):
            data_set_particle = ET.SubElement(collection,"DataSet")
            data_set_particle.attrib['timestep'] = str(self.timesteps[i].currentTime)
            data_set_particle.attrib['group'] = 'Particle centroids'
            data_set_particle.attrib['part'] = '0' 
            data_set_particle.attrib['file'] = str(self.output_particle_file_names[i])
        indent(root_element)         
        tree = ET.ElementTree(root_element)
        tree.write(self.name+'/particles.pvd',encoding="UTF-8",xml_declaration=True)


    def output_vtk_particle_domain_file(self,root_element):
        # All time steps put their valid time steps into a list of times and files
        all_times_files = []
        i = 0
        for t in range(len(self.timesteps)):
           if self.timesteps[t] is not None:
             tsd = { 't': t, 
                     'timestep': self.timesteps[t].currentTime,
                     'file': self.output_particle_domain_file_names[i] }
             i += 1
             all_times_files.append(tsd)

        # Gather all valid timesteps from all tasks
        if self.mpicomm is not None:
          all_times_files = self.mpicomm.gather(all_times_files, root=0)

        if self.mpirank == 0: 
          # Sort all timesteps by time index into new list
          sort_times_files = [None]*len(self.timesteps)
          for rnk in all_times_files:
            for tsd in rnk:
              sort_times_files[tsd['t']] = tsd

          collection = ET.SubElement(root_element,"Collection")
          for i in range(len(self.timesteps)):
              data_set_particle_domain = ET.SubElement(collection,"DataSet")
              data_set_particle_domain.attrib['timestep'] = str(sort_times_files[i]['timestep'])
              data_set_particle_domain.attrib['group'] = 'Particle domains'
              data_set_particle_domain.attrib['part'] = '0'
              data_set_particle_domain.attrib['file'] = str(sort_times_files[i]['file'])
          indent(root_element)         
          tree = ET.ElementTree(root_element)
          tree.write(self.name+'/domains.pvd',encoding="UTF-8",xml_declaration=True)


class Timestep:

    def __init__(self,timestep_xml=None):
        self.dir = ''
        self.datafiles = []
        if timestep_xml is not None:
            Uintah_timestep = self.read_timestep(timestep_xml)
        else:
            return

        meta = Uintah_timestep.find('Meta')
        self.endianness = meta.findtext('endianness')
        self.nBits = int(meta.findtext('nBits'))
        self.numProcs = int(meta.findtext('numProcs'))
        time = Uintah_timestep.find('Time')
        self.timestepNumber = int(time.findtext('timestepNumber'))
        self.currentTime = float(time.findtext('currentTime'))
        self.oldDelt = float(time.findtext('oldDelt'))
        self.grid = Grid(Uintah_timestep.find('Grid'))  
        # All the grid data structures are created, i.e. levels,patches
        # actually read in the l*/p*.xml to read in the actual variable data
        # and store it in each patch's list of variables

        # Read uintah data 
        self.read_datafile(Uintah_timestep)
        self.mat_index = 0 
        self.materials = []
        self.particles_per_material = []
        self.interpolator = ''
        # Extract materials and interpolator from timestep.xml
        self.read_materials(Uintah_timestep)
       
    def read_datafile(self,timestep):
        datafile_iter = timestep.getiterator('Datafile')  #py2.7 use: iter('Datafile')
        for data in datafile_iter:
            da = data.get('href')
            self.datafiles.append(da)

        for files in self.datafiles:
            datafile_xml = self.dir + '/' + files
            datafile_split = datafile_xml.split('/')
            level_id = datafile_split[1]
            tree = ET.parse(datafile_xml)
            elem = tree.getroot()
            for var in elem.getiterator('Variable'):
                v = Variable(var)
                datafile_name = self.dir + '/' + level_id + '/' + v.filename
                v.read_data(datafile_name)
                patch = self.grid.find_patch(v.patch)
                patch.add_variable(v)

    def read_materials(self,timestep):
        # read materials
        tree = ET.parse(self.dir+'/timestep.xml')
        elem = tree.getroot()
        mpm = elem.find('MPM')
        self.interpolator = mpm.find('interpolator').text
        mp_vars = elem.find('MaterialProperties')
        mpm_vars = mp_vars.find('MPM')
        mat_variables_iter = mpm_vars.getiterator('material')
        # determine material numbers
        for mat_vars in mat_variables_iter:
          name = mat_vars.get('name')
          index = mat_vars.get('index')
          self.materials.append(index)
        # Grid uses max index + 1 for all materials
        self.materials.sort()
        self.gmaterials=list(self.materials)
        indexmax=int(max(self.gmaterials))+1
        self.gmaterials.append(str(indexmax))


    def read_timestep(self,timestep_xml):
        self.dir = timestep_xml.split('/')[0]
        try:
            tree = ET.parse(timestep_xml)
            print "   %s" % self.dir 
        except Exception:
            print "Unexpected error opening %s" % timestep_xml 
            return
            
        return tree.getroot()

    def output_vtk_grid(self,filename,patch,mat_index):
        self.grid.output_vtk_grid(filename,patch,mat_index)

    def output_vtk_particle(self,filename,patch,mat_index):
        self.grid.output_vtk_particle(filename,patch,mat_index)

    def output_vtk_particle_domain(self,filename,patch,mat_index):
        self.grid.output_vtk_particle_domain(filename,patch,mat_index,self.interpolator)

    def output_vtk_pvti_grid(self,filename,subfiles):
        self.grid.output_vtk_pvti_grid(filename,subfiles)

    def output_vtk_pvtu_particle(self,filename,subfiles):
        self.grid.output_vtk_pvtu_particle(filename,subfiles)

    def output_vtk_pvtu_domain(self,filename,subfiles):
        self.grid.output_vtk_pvtu_domain(filename,subfiles)


class Grid:
    def __init__(self,grid):
        self.numLevels = int(grid.findtext('numLevels'))
        self.levels = []
        level_iter = grid.getiterator('Level')
        self.read_level(level_iter)

    def print_grid(self):
        print "Grid number of levels = %s" % self.numLevels

    def read_level(self,level_iter):
        for level in level_iter:
            self.add_level(Level(level))

    def add_level(self,level):
        self.levels.append(level)

        
    def get_levels(self):
        return self.levels

    def get_grid_extents(self):
        extent = self.get_extent()
        print "Grid extent = %s %s" % (extent[0],extent[1])
        return extent

    def get_extent(self):
        ex = self.levels[0].get_extent()
        for level in self.levels:
            extent = level.get_extent()
            # print "Level %s extent = %s" % (level.id,extent)
            ex_lo = [min(ex[0][0],extent[0][0]),min(ex[0][1],extent[0][1]), \
                         min(ex[0][2],extent[0][2])]
            ex_hi =[max(ex[1][0],extent[1][0]), max(ex[1][1],extent[1][1]), \
                        max(ex[1][2],extent[1][2])]
            # print "Current ex_lo = %s" % ex_lo
            # print "Current ex_hi = %s" % ex_hi
            ex = (ex_lo,ex_hi)

        return ex


    def find_patch(self,id):
        for level in self.levels:
            for patch in level.patches:
                if id == patch.id:
                    return patch
            
        return None

    def vtk_grid(self,patch,endianness=None,nbits=None):
        # Not using patch info yet?
        extent = self.get_grid_extents()
        lo = extent[0]
        hi = extent[1]
        string_extent = repr(lo[0]) + ' ' + repr(hi[0]) \
            + ' ' + repr(lo[1]) + ' ' + repr(hi[1]) \
            + ' ' + repr(lo[2]) + ' ' + repr(hi[2]) 
        # rectilinear_elem = create_vtkfile_element(vtk_type='RectilinearGrid')
        # grid_elem = ET.Element('RectilinearGrid')
        rectilinear_elem = create_vtkfile_element(vtk_type='ImageData')
        grid_elem = ET.Element('ImageData')
        grid_elem.attrib['WholeExtent'] = string_extent
        for level in self.levels:
            level.vtk_element(patch,grid_elem)
        rectilinear_elem.append(grid_elem)
        indent(rectilinear_elem)
        #ET.dump(rectilinear_elem)
        return rectilinear_elem


    def output_vtk_grid(self,filename,patch,mat_index):
        for level in self.levels:
            elem = level.vtk_grid(patch,mat_index)
            elem_tree = ET.ElementTree(elem)
            elem_tree.write(filename,encoding="UTF-8",xml_declaration=True)

    def output_vtk_particle(self,filename,patch,mat_index):
        for level in self.levels:
            elem = level.vtk_particle(patch,mat_index)
            elem_tree = ET.ElementTree(elem)
            elem_tree.write(filename,encoding="UTF-8",xml_declaration=True)


    def output_vtk_particle_domain(self,filename,patch,mat_index,interpolator):
        for level in self.levels:
            elem = level.vtk_particle_domain(patch,mat_index,interpolator)
            elem_tree = ET.ElementTree(elem)
            #ET.dump(elem)
            elem_tree.write(filename,encoding="UTF-8",xml_declaration=True)

    def output_vtk_pvti_grid(self,filename,subfiles):
        for level in self.levels:
            elem = level.vtk_pvti_grid(subfiles)
            elem_tree = ET.ElementTree(elem)
            elem_tree.write(filename,encoding="UTF-8",xml_declaration=True)
    
    def output_vtk_pvtu_particle(self,filename,subfiles):
        for level in self.levels:
            elem = level.vtk_pvtu_particle(subfiles)
            elem_tree = ET.ElementTree(elem)
            elem_tree.write(filename,encoding="UTF-8",xml_declaration=True)
    
    def output_vtk_pvtu_domain(self,filename,subfiles):
        for level in self.levels:
            elem = level.vtk_pvtu_domain(subfiles)
            elem_tree = ET.ElementTree(elem)
            elem_tree.write(filename,encoding="UTF-8",xml_declaration=True)


class Level:
    def __init__(self,level):
        self.numPatches = int(level.findtext('numPatches'))
        self.totalCells = int(level.findtext('totalCells'))
        self.extraCells = IntVector(level.findtext('extraCells'))
        self.anchor = Vector(level.findtext('anchor'))
        self.id = int(level.findtext('id'))
        self.cellspacing = Vector(level.findtext('cellspacing'))
        self.patches = []
        self.varlist = []
        patch_iter = level.getiterator('Patch')
        self.read_patch(patch_iter)
        for p in self.patches:
            p.find_neighbors(self.patches)


    def print_level(self):
        print "Level number of patches = %s" % self.numPatches
        print "Level number of total cells  = %s" % self.totalCells
        print "Level extra cells = %s" % self.extraCells
        print "Level anchor = %s" % self.anchor
        print "Level id = %s" % self.id
        print "Level cell spacing = %s" % self.cellspacing


    def read_patch(self,patch_iter):
        for patch in patch_iter:
            self.add_patch(Patch(patch))

        
    def add_patch(self,patch):
        self.patches.append(patch)

    
    def get_patches(self):
        return self.patches


    def get_extent(self):
        ex = self.patches[0].get_extent()
        for patch in self.patches:
            extent = patch.get_extent()
            # print "Patch %s extent = %s" % (patch.id,extent)
            ex_lo = [min(ex[0][0],extent[0][0]),min(ex[0][1],extent[0][1]), \
                         min(ex[0][2],extent[0][2])]
            ex_hi =[max(ex[1][0],extent[1][0]), max(ex[1][1],extent[1][1]), \
                        max(ex[1][2],extent[1][2])]
            # print "Current ex_lo = %s" % ex_lo
            # print "Current ex_hi = %s" % ex_hi
            ex = (ex_lo,ex_hi)

        return ex


    def vtk_grid(self,patch,mat_index):
        extent = self.get_extent()
        #extent = patch.get_extent()
        lo = extent[0]
        hi = extent[1]
        string_extent = repr(lo[0]) + ' ' + repr(hi[0]) \
            + ' ' + repr(lo[1]) + ' ' + repr(hi[1]) \
            + ' ' + repr(lo[2]) + ' ' + repr(hi[2]) 
        vtkfile_elem = create_vtkfile_element(vtk_type='ImageData')
        image_data_elem = ET.Element('ImageData')
        image_data_elem.attrib['WholeExtent'] = string_extent
        image_data_elem.attrib['Origin'] = repr(self.anchor[0]) + ' ' +\
            repr(self.anchor[1]) + ' ' + repr(self.anchor[2])
        image_data_elem.attrib['Spacing'] = repr(self.cellspacing[0]) + ' ' +\
            repr(self.cellspacing[1]) + ' ' + repr(self.cellspacing[2])
        #for patch in self.patches:
        #    elem = patch.vtk_particle_domain(unstructured_data_elem,mat_index,interpolator)
        #    unstructured_data_elem.append(elem)
        elem = patch.vtk_grid(image_data_elem,mat_index)
        image_data_elem.append(elem)
        vtkfile_elem.append(image_data_elem)
        #vtkfile_elem.append(elem)
        indent(vtkfile_elem)
        #ET.dump(vtkfile_elem)
        return vtkfile_elem

    def vtk_particle(self,patch,mat_index):
        vtkfile_elem = create_vtkfile_element(vtk_type='UnstructuredGrid')
        unstructured_data_elem = ET.Element('UnstructuredGrid')
        elem = patch.vtk_particle(unstructured_data_elem,mat_index)
        unstructured_data_elem.append(elem)
        vtkfile_elem.append(unstructured_data_elem)
        indent(vtkfile_elem)
        return vtkfile_elem

    def vtk_particle_domain(self,patch,mat_index,interpolator):
        vtkfile_elem = create_vtkfile_element(vtk_type='UnstructuredGrid')
        unstructured_data_elem = ET.Element('UnstructuredGrid')
        elem = patch.vtk_particle_domain(unstructured_data_elem,mat_index,interpolator)
        unstructured_data_elem.append(elem)
        vtkfile_elem.append(unstructured_data_elem)
        indent(vtkfile_elem)
        return vtkfile_elem

    def vtk_pvti_grid(self,subfiles):
        #self.print_level()
        # need to fix grid data (Check ghost levels, what to do for all materials, etc.)
        extent = self.get_extent()
        lo = extent[0]
        hi = extent[1]
        string_extent = repr(lo[0]) + ' ' + repr(hi[0]) \
            + ' ' + repr(lo[1]) + ' ' + repr(hi[1]) \
            + ' ' + repr(lo[2]) + ' ' + repr(hi[2]) 
        vtkfile_elem = create_vtkfile_element(vtk_type='PImageData')
        pimage_data_elem = ET.Element('PImageData')
        pimage_data_elem.attrib['WholeExtent'] = string_extent
        pimage_data_elem.attrib['GhostLevel'] = str(self.extraCells[0])
        pimage_data_elem.attrib['Origin'] = repr(self.anchor[0]) + ' ' +\
            repr(self.anchor[1]) + ' ' + repr(self.anchor[2])
        pimage_data_elem.attrib['Spacing'] = repr(self.cellspacing[0]) + ' ' +\
            repr(self.cellspacing[1]) + ' ' + repr(self.cellspacing[2])
        # may want all materials on grid later
        for patch in self.patches:  
            elem = patch.vtk_pvti_grid(pimage_data_elem,subfiles)
        # One list of grid variables for all patches
        pimage_data_elem = self.read_vti_subfiles(elem,subfiles)
        indent(pimage_data_elem)
        #vtkfile_elem.extend(elem)
        vtkfile_elem.append(pimage_data_elem)
        indent(vtkfile_elem)
        #ET.dump(vtkfile_elem)
        return vtkfile_elem
    

    def vtk_pvtu_particle(self,subfiles):
        vtkfile_elem = create_vtkfile_element(vtk_type='PUnstructuredGrid')
        punstructured_data_elem = ET.Element('PUnstructuredGrid')
        for patch in self.patches:
            elem = patch.vtk_pvtu_particle(punstructured_data_elem,subfiles)
            punstructured_data_elem = self.read_vtu_subfiles(elem,subfiles)
        indent(punstructured_data_elem)
        vtkfile_elem.append(punstructured_data_elem)
        indent(vtkfile_elem)
        return vtkfile_elem

    def vtk_pvtu_domain(self,subfiles):
        vtkfile_elem = create_vtkfile_element(vtk_type='PUnstructuredGrid')
        punstructured_data_elem = ET.Element('PUnstructuredGrid')
        #rbl:  Use for grid spacing.
        #      print "(level) grid spacing=",self.cellspacing
        for patch in self.patches:
            elem = patch.vtk_pvtu_domain(punstructured_data_elem,subfiles)
            punstructured_data_elem = self.read_vtu_subfiles(elem,subfiles)
        indent(punstructured_data_elem)
        vtkfile_elem.append(punstructured_data_elem)
        indent(vtkfile_elem)
        #ET.dump(vtkfile_elem)
        return vtkfile_elem

    def read_vti_subfiles(self,vtkfile_elem,subfiles):
        ''' 
        Reads created material summary grid subfiles and creates master lists of DataArrays.
        Assumes that all patches and materials have the same number and list of variables.
        '''
        #points_attrib = []
        #cells_attrib = []
        celldata_attrib = []
        pointdata_attrib = []
        # Don't print some elements if ImageData file
        #image = [s for s in subfiles if "vti" in s]
        #if not image:
        #  points = ET.SubElement(vtkfile_elem,'PPoints')        
        #  cells = ET.SubElement(vtkfile_elem,'PCells')        
        pointdata = ET.SubElement(vtkfile_elem,'PPointData')        
        celldata = ET.SubElement(vtkfile_elem,'PCellData')        
        for sub in subfiles:
          tree = ET.ElementTree(file=sub)
          root = tree.getroot()
          for elem in tree.iter():
              if elem.tag == 'CellData':
                for subelem in elem:
                  if subelem.attrib not in celldata_attrib:
                    subelem.tag = 'P'+subelem.tag
                    subelem.text=''
                    celldata_attrib.append(subelem.attrib)
                    celldata.append(subelem)
              elif elem.tag == 'PointData':
                for subelem in elem:
                  if subelem.attrib not in pointdata_attrib:
                    subelem.tag = 'P'+subelem.tag
                    subelem.text=''
                    pointdata_attrib.append(subelem.attrib)
                    pointdata.append(subelem)
              elif elem.tag == 'Points':
                for subelem in elem:
                  if subelem.attrib not in points_attrib:
                    subelem.tag = 'P'+subelem.tag
                    subelem.text=''
                    points_attrib.append(subelem.attrib)
                    points.append(subelem)
              elif elem.tag == 'Cells':
                for subelem in elem:
                  if subelem.attrib not in cells_attrib:
                    subelem.tag = 'P'+subelem.tag
                    subelem.text=''
                    cells_attrib.append(subelem.attrib)
                    cells.append(subelem)
        indent(vtkfile_elem)
        return vtkfile_elem 


    def read_vtu_subfiles(self,vtkfile_elem,subfiles):
        ''' 
        Reads created material particle subfiles and creates master lists of DataArrays.
        Assumes that all materials and patches have the same number and list of variables.
        '''
        points_attrib = []
        cells_attrib = []
        celldata_attrib = []
        pointdata_attrib = []
        points = ET.SubElement(vtkfile_elem,'PPoints')        
        cells = ET.SubElement(vtkfile_elem,'PCells')        
        pointdata = ET.SubElement(vtkfile_elem,'PPointData')        
        celldata = ET.SubElement(vtkfile_elem,'PCellData')        
        for sub in subfiles:
          tree = ET.ElementTree(file=sub)
          root = tree.getroot()
          for elem in tree.iter():
              if elem.tag == 'CellData':
                for subelem in elem:
                  if subelem.attrib not in celldata_attrib:
                    subelem.tag = 'P'+subelem.tag
                    subelem.text=''
                    celldata_attrib.append(subelem.attrib)
                    celldata.append(subelem)
              elif elem.tag == 'PointData':
                for subelem in elem:
                  if subelem.attrib not in pointdata_attrib:
                    subelem.tag = 'P'+subelem.tag
                    subelem.text=''
                    pointdata_attrib.append(subelem.attrib)
                    pointdata.append(subelem)
              elif elem.tag == 'Points':
                for subelem in elem:
                  if subelem.attrib not in points_attrib:
                    subelem.tag = 'P'+subelem.tag
                    subelem.text=''
                    points_attrib.append(subelem.attrib)
                    points.append(subelem)
              elif elem.tag == 'Cells':
                for subelem in elem:
                  if subelem.attrib not in cells_attrib:
                    subelem.tag = 'P'+subelem.tag
                    subelem.text=''
                    cells_attrib.append(subelem.attrib)
                    cells.append(subelem)
        indent(vtkfile_elem)
        return vtkfile_elem 


class Patch:
    def __init__(self,patch):
        self.id = int(patch.findtext('id'))
        self.proc = int(patch.findtext('proc'))
        self.lowIndex = IntVector(patch.findtext('lowIndex'))
        self.highIndex = IntVector(patch.findtext('highIndex'))
        self.interiorLowIndex = IntVector(patch.findtext('interiorLowIndex'))
        self.interiorHighIndex = IntVector(patch.findtext('interiorHighIndex'))
        self.nnodes = int(patch.findtext('nnodes'))
        self.lower = Vector(patch.findtext('lower'))
        self.upper = Vector(patch.findtext('upper'))
        self.totalCells = int(patch.findtext('totalCells'))
        self.variables = []
        self.varlist = []
        self.plus_neighbor = [0,0,0]


    def print_patch(self):
        print "Patch id = %s" % self.id
        print "Patch proc = %s" % self.proc
        print "Patch lowIndex = %s" % self.lowIndex
        print "Patch highIndex = %s" % self.highIndex
        print "Patch interiorLowIndex = %s" % self.interiorLowIndex
        print "Patch interiorHighIndex = %s" % self.interiorHighIndex
        print "Patch nnodes = %s" % self.nnodes
        print "Patch lower = %s" % self.lower
        print "Patch upper = %s" % self.upper
        print "Patch totalCells = %s" % self.totalCells


    def add_variable(self,variable):
        self.variables.append(variable)
    
     
    def get_variables(self):
        return self.variables

    
    def generate_grid(self):
        pass


    def find_neighbors(self,neighbors):
        for n in neighbors:
            if self.id == n.id:
                continue
            else:
                if self.highIndex[0] == n.lowIndex[0]:
                    self.plus_neighbor[0] = 1
                if self.highIndex[1] == n.lowIndex[1]:
                    self.plus_neighbor[1] = 1
                if self.highIndex[2] == n.lowIndex[2]:
                    self.plus_neighbor[2] = 1
                
    def get_extent(self):
        return (self.lowIndex, self.highIndex)

    def get_material_list(self):
        material = [] 
        for v in self.variables:
            if v.index not in material:
              material.append(v.index) 
        material.sort()
        return material

    def get_num_materials(self):
        num_mat = 0
        for v in self.variables:
            if v.type == 'ParticleVariable<Point>':
                num_mat += 1
        return num_mat


    def get_num_particles(self,mat_index):
        for v in self.variables:
            if v.type == 'ParticleVariable<Point>':
                if v.index == int(mat_index):
                    #print 'mat=',mat_index,' v.numParticles=',v.numParticles
                    return v.numParticles 
 
    def get_interpolator(self):
        return self.interpolator

    def get_particle_list(self, mat_index):
        particles = [] 
        for v in self.variables:
            if (v.type == 'ParticleVariable<Point>' and v.index == mat_index):
              if v.type not in particles:
                particles.append(v.index) 
        #particles.sort()
        return particles 

    def file_patch_numbers(self,subfiles):
        # break material files into separate patch files
        newfiles=[]
        for s in subfiles:
          tmp = s.split('.')
          ext = '.'.join(tmp[-1:])
          patchid= tmp[-2].split('_p')[-1]
          if self.id == int(patchid):
            return s

    def subfile_patch_number(self,subfile):
        # return corresponding patch subfile
        newfile=[]
        tmp = subfile.split('.')
        ext = '.'.join(tmp[-1:])
        patchid= tmp[-2].split('_p')[-1]
        if self.id == int(patchid):
            return True 


    def vtk_pvti_grid(self,root_elem,subfiles):
        extent = self.get_extent()
        extent = self.get_extent()
        lo = extent[0]
        hi = extent[1]
        string_extent = repr(lo[0]) + ' ' + repr(hi[0]) \
          + ' ' + repr(lo[1]) + ' ' + repr(hi[1]) \
          + ' ' + repr(lo[2]) + ' ' + repr(hi[2]) 
        # PPointData
        #pointdata_elem = ET.SubElement(root_elem,'PPointData')
        #pointdata_elem.attrib['Scalars'] = 'Material'
        #pointdata_array = ET.SubElement(pointdata_elem,'PDataArray')
        #pointdata_array.attrib['type'] = 'Int32'
        #pointdata_array.attrib['Name'] = 'Material' 
        #pointdata_array.attrib['NumberOfComponents'] = '1' 
        # PCellData 
        #celldata_elem = ET.SubElement(root_elem,'PCellData')
        # Loop over patches and materials 
        patchid = self.id
        for sub in subfiles:
          if self.subfile_patch_number(sub):
            piece_elem = ET.Element('Piece')
            piece_elem.attrib['Source']=sub
            piece_elem.attrib['Extent']=string_extent
            root_elem.append(piece_elem)
            #ET.dump(root_elem) 
            return root_elem


    def vtk_pvtu_particle(self,root_elem,subfiles):
        patch_elem = ET.SubElement(root_elem,'PUnstructuredGrid')
        # rbl: replace with actual ghost level from uintah if needed for particles?
        patch_elem.attrib['GhostLevel'] = '0'
        # PPoints
        points_elem = ET.Element('PPoints')
        pointsdata_array = ET.SubElement(points_elem,'PDataArray')
        pointsdata_array.attrib['type'] = 'Int32'
        pointsdata_array.attrib['Name'] = 'Material' 
        pointsdata_array.attrib['NumberOfComponents'] = '1' 
        # PCells
        cells_elem = ET.Element('PCells')
        #cells_elem.attrib['type'] = 'Float64'
        cells_elem.attrib['type'] = 'Float32'
        cells_elem.attrib['Name'] = 'Cells' 
        cells_elem.attrib['NumberOfComponents'] = '3' 
        cells_elem.attrib['format'] = 'ascii'
        # Cell connectivity
        cell_connectivity = ET.SubElement(cells_elem,'PDataArray')
        cell_connectivity.attrib['type'] = 'Int32'
        cell_connectivity.attrib['Name'] = 'connectivity' 
        cell_connectivity.attrib['format'] = 'ascii'
        # Cell offsets 
        cell_offsets = ET.SubElement(cells_elem,'PDataArray')
        cell_offsets.attrib['type'] = 'Int32'
        cell_offsets.attrib['Name'] = 'offsets' 
        cell_offsets.attrib['NumberOfComponents'] = '1' 
        cell_offsets.attrib['format'] = 'ascii'
        #cell_offsets.text=''
        # Cell types 
        cell_types = ET.SubElement(cells_elem,'PDataArray')
        cell_types.attrib['type'] = 'UInt8'
        cell_types.attrib['Name'] = 'types' 
        cell_types.attrib['NumberOfComponents'] = '1' 
        cell_types.attrib['format'] = 'ascii'
        # PPointData 
        pointdata_elem = ET.Element('PPointData')
        # PCellData 
        celldata_elem = ET.Element('PCellData')
        # loop over materials and patches
        for sub in subfiles:
          piece_elem = ET.Element('Piece')
          piece_elem.attrib['Source']=sub
          patch_elem.append(piece_elem)
        indent(patch_elem)
        #print 'patch_elem=',ET.dump(patch_elem)
        return patch_elem



    def vtk_pvtu_domain(self,root_elem,subfiles):
        patch_elem = ET.SubElement(root_elem,'PUnstructuredGrid')
        # rbl: replace with actual ghost level from uintah?
        patch_elem.attrib['GhostLevel'] = '0'
        # PPoints
        points_elem = ET.Element('PPoints')
        # PPointData 
        pointdata_elem = ET.Element('PPointData')
        # PCells
        cells_elem = ET.Element('PCells')
        #cells_elem.attrib['type'] = 'Float64'
        cells_elem.attrib['type'] = 'Float32'
        cells_elem.attrib['Name'] = 'Cells' 
        cells_elem.attrib['NumberOfComponents'] = '3' 
        cells_elem.attrib['format'] = 'ascii'
        # Cell connectivity
        cell_connectivity = ET.SubElement(cells_elem,'PDataArray')
        cell_connectivity.attrib['type'] = 'Int32'
        cell_connectivity.attrib['Name'] = 'connectivity' 
        cell_connectivity.attrib['format'] = 'ascii'
        # Cell offsets 
        cell_offsets = ET.SubElement(cells_elem,'PDataArray')
        cell_offsets.attrib['type'] = 'Int32'
        cell_offsets.attrib['Name'] = 'offsets' 
        cell_offsets.attrib['NumberOfComponents'] = '1' 
        cell_offsets.attrib['format'] = 'ascii'
        cell_offsets.text=''
        # Cell types 
        cell_types = ET.SubElement(cells_elem,'PDataArray')
        cell_types.attrib['type'] = 'UInt8'
        cell_types.attrib['Name'] = 'types' 
        cell_types.attrib['NumberOfComponents'] = '1' 
        cell_types.attrib['format'] = 'ascii'
        # PCellData 
        celldata_elem = ET.Element('PCellData')
        # loop over patches and materials
        for sub in subfiles:
          piece_elem = ET.Element('Piece')
          piece_elem.attrib['Source']=sub
          patch_elem.append(piece_elem)
        indent(patch_elem)
        #print 'patch_elem=',ET.dump(patch_elem)
        return patch_elem


    def vtk_grid(self,root_elem,mat_index):
        num_nodes = self.nnodes
        #num_elements = self.nnodes
        #print 'self.nnodes=',self.nnodes
        #print 'self.totalCells=',self.totalCells
        extent = self.get_extent()
        lo = extent[0]
        hi = extent[1]
        hi = map(operator.sub,hi,self.plus_neighbor)
        # subtract off the plus_neighbor values from the hi extent
        # hi = extent[1] - plus_neighbor
        string_extent = repr(lo[0]) + ' ' + repr(hi[0]) \
            + ' ' + repr(lo[1]) + ' ' + repr(hi[1]) \
            + ' ' + repr(lo[2]) + ' ' + repr(hi[2]) 
        patch_elem = ET.Element('Piece')
        patch_elem.attrib['Extent'] = string_extent

        pointdata_elem = ET.Element('PointData')
        #pointdata_array = ET.SubElement(pointdata_elem,'DataArray')
        # rbl: changed to just putting all materials toghether on grid
        #pointdata_array.attrib['Name'] = 'Material' 
        #pointdata_array.attrib['type'] = 'Int32'
        #pointdata_array.attrib['NumberOfComponents'] = '1' 

        celldata_elem = ET.Element('CellData')
        for variable in self.variables:
            variable_type = variable.type.split('<')[0]
            mat = variable.index
            elem = variable.vtk_element()
            #print 'elem=',elem.text 
            # CC 1 value, NC 8 (Uintah API says 4?), FC 6
            #if variable_type == 'CCVariable' and mat == int(mat_index):
            #if variable_type == 'FCVariable' and mat == int(mat_index):
            #    celldata_elem.append(elem)
            if variable_type == 'NCVariable' and mat == int(mat_index):
                pointdata_elem.append(elem)
        # Convert to text 
        pointdata_text = ''
        celltypes_text = ''
        celloffsets_text = ''
        cellconnectivity_text = ''
        mat_text = ''
        dim = 3
        # perhaps need a new reference for grid nodes or elements
        #for g in range(0,num_nodes):
        #    # append offset (from end of first element)
        #    celloffsets_text += str(dim*(g+1))
        #    # append type 
        #    #celltypes_text += vtk_cell_type 
        #    # append material 
        #    mat_text += mat_index+' ' 
        # Store results
        #pointdata_array.text = mat_text
        #points_elem.text=pointdata_text
        #mat_text += mat_index+' ' 
        patch_elem.append(pointdata_elem)
        patch_elem.append(celldata_elem)
        indent(patch_elem)
        #ET.dump(patch_elem)
        return patch_elem
 


    def vtk_particle(self,root_elem,mat_index):
          patch_elem = ET.Element('Piece')
          num_particles = self.get_num_particles(mat_index)
          #   1     Point           VTK_VERTEX      Particle centroid
          vtk_cell_type = '1 '                                                  # Vertex only (Particle centroid) 
          patch_elem.attrib['NumberOfPoints'] = repr(num_particles)
          patch_elem.attrib['NumberOfCells'] =  '0' 
          # Cells
          cells_elem = ET.Element('Cells')
          # Cell connectivity
          cell_connectivity = ET.SubElement(cells_elem,'DataArray')
          cell_connectivity.attrib['type'] = 'Int32'
          cell_connectivity.attrib['Name'] = 'connectivity' 
          cell_connectivity.attrib['format'] = 'ascii'
          # Cell offsets 
          cell_offsets = ET.SubElement(cells_elem,'DataArray')
          cell_offsets.attrib['type'] = 'Int32'
          cell_offsets.attrib['Name'] = 'offsets' 
          cell_offsets.attrib['format'] = 'ascii'
          # Cell types 
          cell_types = ET.SubElement(cells_elem,'DataArray')
          cell_types.attrib['type'] = 'UInt8'
          cell_types.attrib['Name'] = 'types' 
          cell_types.attrib['format'] = 'ascii'
          # Cell data
          celldata_elem = ET.Element('CellData')
          # Points (extracted from variable loop below)
          points_elem = ET.Element('Points')
          # Point data
          pointdata_elem = ET.Element('PointData')
          pointdata_elem.attrib['Scalars'] = 'Material' 
          pointdata_array = ET.SubElement(pointdata_elem,'DataArray')
          pointdata_array.attrib['Name'] = 'Material' 
          pointdata_array.attrib['type'] = 'Int32'
          pointdata_array.attrib['NumberOfComponents'] = '1' 
          dim = 3                                                               # dimension 
          getCentroid = True                                                    # flags for getting particle lists
          centroid=[] 
          for variable in self.variables:
            variable_type = variable.type.split('<')[0]
            mat = variable.index
            elem = variable.vtk_element()
            if variable_type == 'ParticleVariable' and mat == int(mat_index):
                if variable.data_type == 'Point':
                    points_elem.append(elem)
                    if getCentroid:
                      getCentroid=False
                      centroid =list(chunks(elem.text.split(),vector))
                else:
                    pointdata_elem.append(elem)
          # Convert to text 
          pointdata_text = ''
          celltypes_text = ''
          celloffsets_text = ''
          cellconnectivity_text = ''
          mat_text = ''
          for p in range(0,num_particles):
            # append offset (from end of first element)
            celloffsets_text += str(dim*(p+1))
            # append type 
            celltypes_text += vtk_cell_type 
            # append material 
            mat_text += mat_index+' ' 
          # Store results
          pointdata_array.text = mat_text
          points_elem.text=pointdata_text
          cell_connectivity.text=cellconnectivity_text
          cell_offsets.text=celloffsets_text
          cell_types.text=celltypes_text
          # Append patch
          patch_elem.append(points_elem)
          patch_elem.append(cells_elem)
          patch_elem.append(pointdata_elem)
          patch_elem.append(celldata_elem)
          indent(patch_elem)
          return patch_elem



    def vtk_particle_domain(self,root_elem,mat_index,interpolator):
          patch_elem = ET.Element('Piece')
          num_particles = self.get_num_particles(mat_index)
          #   1     Point           VTK_VERTEX      Particle centroid
          #   10    Tetrahedron     VTK_TETRA       Cpti particle domain
          #   11    Parallelpiped   VTK_VOXEL       Mpm particle domain
          #   12    Hexahedron      VTK_HEXAHEDRON  Cpdi/Cbdi particle domain
          if interpolator == 'cpti':
            cpp = 4                                                             # CPTI 4 element corners 
            vtk_cell_type = '10 '                                               # Tetrahedron
          else:
            cpp = 8                                                             # CPDI,CBDI, GIMP 8 element corners 
            vtk_cell_type = '12 '                                               # Hexahedron
          num_points = num_particles*cpp                                        # total number of corners
          patch_elem.attrib['NumberOfPoints'] = repr(num_points)
          patch_elem.attrib['NumberOfCells'] =  repr(num_particles)
          # Cells
          cells_elem = ET.Element('Cells')
          # Cell connectivity
          cell_connectivity = ET.SubElement(cells_elem,'DataArray')
          cell_connectivity.attrib['type'] = 'Int32'
          cell_connectivity.attrib['Name'] = 'connectivity'
          cell_connectivity.attrib['format'] = 'ascii'
          # Cell offsets 
          cell_offsets = ET.SubElement(cells_elem,'DataArray')
          cell_offsets.attrib['type'] = 'Int32'
          cell_offsets.attrib['Name'] = 'offsets'
          cell_offsets.attrib['format'] = 'ascii'
          # Cell types 
          cell_types = ET.SubElement(cells_elem,'DataArray')
          cell_types.attrib['type'] = 'UInt8'
          cell_types.attrib['Name'] = 'types'
          cell_types.attrib['format'] = 'ascii'
          # Cell data
          celldata_elem = ET.Element('CellData')
          celldata_elem.attrib['Scalars'] = 'Material'
          celldata_array = ET.SubElement(celldata_elem,'DataArray')
          celldata_array.attrib['Name'] = 'Material'
          celldata_array.attrib['type'] = 'Int32'
          celldata_array.attrib['NumberOfComponents'] = '1'
          # Points
          points_elem = ET.Element('Points')
          points_data_elem = ET.SubElement(points_elem,'DataArray')
          points_data_elem.attrib['type'] = 'Float32'
          points_data_elem.attrib['Name'] = 'Particle Domains'
          points_data_elem.attrib['NumberOfComponents'] = '3'
          # PointData
          pointdata_elem = ET.Element('PointData')
          dim = 3                                                               # dimension 
          gotCentroid = False                                                   # flags for getting particle lists
          gotPscalefactor = False
          gotPsize = False
          gotDefgrad = False
          xpc=np.zeros((num_particles,cpp,dim),float)                           # particle domain corners
          r1=np.zeros((num_particles,dim),float)                                # particle Rvectors 
          r2=np.zeros((num_particles,dim),float)                                # particle Rvectors 
          r3=np.zeros((num_particles,dim),float)                                # particle Rvectors 
          centroid=[]
          pSize=[]
          pscalefactor=[]
          defgrad=[]
          for variable in self.variables:
            variable_type = variable.type.split('<')[0]
            mat = variable.index
            elem = variable.vtk_element()
            # perhaps map or some quicker breakout of materials?
            if variable_type == 'ParticleVariable' and mat == int(mat_index):
                if variable.data_type == 'Point':
                    #points_elem.append(elem)
                    if not gotCentroid:
                      gotCentroid=True
                      centroid =list(chunks(elem.text.split(),vector))
                elif variable.data_type == 'Matrix3':
                    celldata_elem.append(elem)
                    if variable.name == 'p.scalefactor' and not gotPscalefactor:
                      gotPscalefactor=True
                      pscalefactor = list(chunks(elem.text.split(),tensor))
                    elif variable.name == 'p.size' and not gotPsize:
                      gotPsize=True
                      pSize = list(chunks(elem.text.split(),tensor))
                    elif variable.name == 'p.deformationMeasure' and not gotDefgrad:
                      gotDefgrad=True
                      defgrad =list(chunks(elem.text.split(),tensor))
                else:
                    celldata_elem.append(elem)
          # Extract particle domain corners and connectivity
          pointdata_text = ''
          celltypes_text = ''
          celloffsets_text = ''
          cellconnectivity_text = ''
          mat_text = ''
          for p in range(0,num_particles):
            if gotPscalefactor:
              r1[p] = (pscalefactor[p][0],pscalefactor[p][3],pscalefactor[p][6])
              r2[p] = (pscalefactor[p][1],pscalefactor[p][4],pscalefactor[p][7])
              r3[p] = (pscalefactor[p][2],pscalefactor[p][5],pscalefactor[p][8])
              # append particle domain corners (perhaps need to have identity in defgrad and pSize?)
              xpc[p] = calcDomain(centroid[p], r1[p], r2[p], r3[p],interpolator)
              #print "particle domains from p.scalefactor"
            elif args.do_size and (gotPsize or gotDefgrad):
              if gotPsize:
                r1[p] = (pSize[p][0],pSize[p][3],pSize[p][6])
                r2[p] = (pSize[p][1],pSize[p][4],pSize[p][7])
                r3[p] = (pSize[p][2],pSize[p][5],pSize[p][8])
              else:
                r1[p] = (1.0,0.0,0.0)
                r2[p] = (0.0,1.0,0.0)
                r3[p] = (0.0,0.0,1.0)
              if not gotDefgrad:
                defgrad[p] = (1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0)
              # need to get grid spacing from level.cellspacing
              xpc[p] = oldCalcDomain(centroid[p], defgrad[p], r1[p], r2[p], r3[p],interpolator)
              print "particle domains from p.size, p.deformationMeasure and grid spacing"
            else:
              print "missing required files for particle domains"

            for c in range(0,cpp):
              corners= ' '.join(map(str,xpc[p][c]))
              pointdata_text += corners+' '
              cellconnectivity_text += str(cpp*p+c)+' '
            pointdata_text += '\n'
            # append connectivity 
            cellconnectivity_text += '\n'
            # append offset (from end of first element)
            celloffsets_text += str(cpp*(p+1))+' '
            # append type 
            celltypes_text += vtk_cell_type
            # append material 
            mat_text += mat_index+' '
          # String results
          celldata_array.text = mat_text
          points_data_elem.text=pointdata_text
          cell_connectivity.text=cellconnectivity_text
          cell_offsets.text=celloffsets_text
          cell_types.text=celltypes_text
          # Append patch
          patch_elem.append(points_elem)
          patch_elem.append(cells_elem)
          patch_elem.append(celldata_elem)
          patch_elem.append(pointdata_elem)
          indent(patch_elem)
          #ET.dump(celldata_array)
          return patch_elem



class Variable:
    def __init__(self,variable):
        self.type = variable.get('type')
        self.data_type = find_type(self.type)
        self.name = variable.findtext('variable')
        self.index = int(variable.findtext('index'))
        self.patch = int(variable.findtext('patch'))
        self.start = int(variable.findtext('start'))
        numParticles = variable.findtext('numParticles')
        if numParticles is not None:
          numParticles = int(numParticles)
        self.numParticles = numParticles
        self.end = int(variable.findtext('end'))
        self.filename = variable.findtext('filename')


    def print_variable(self):
        print "Variable type = %s" % self.type
        print "Variable data type = %s" % self.data_type
        print "Variable name = %s" % self.name
        print "Variable index = %s" % self.index
        print "Variable patch = %s" % self.patch
        print "Variable start = %s" % self.start
        print "Variable numParticles = %s" % self.numParticles
        print "Variable end = %s" % self.end
        print "Variable filename = %s" % self.filename

    def read_data(self,filename):
        fd = open(filename,'rb')
        fd.seek(self.start)
        #        print 'location of file descriptor after start = %s' % fd.tell()
        bytes_per_data = 8
        num_data = (self.end - self.start)/bytes_per_data
        #        print 'Number of data items to load = %d' % num_data
        data_per_node = 1
        if self.data_type == 'Vector' or self.data_type == 'Point':
            data_per_node = 3
        if self.data_type == 'Matrix3':
            data_per_node = 9
        # p.particleID and other large integers shouldn't be loaded as floats 
        if self.type == 'ParticleVariable<long64>':
          self.data = np.fromfile(fd,int,num_data)
          #print str(self.name)+' '+str(self.type)+'='+str(self.data)
        else:
          self.data = np.fromfile(fd,float,num_data)
        fd.close()
        # print self.data.shape
        # print num_data/data_per_node,data_per_node
        self.data = self.data.reshape(num_data/data_per_node,data_per_node)
        
    def get_data(self):
        return self.data


    def vtk_element(self):
        variable_type = self.type.split('<')[0]
        number_type = (self.type.split('<')[1]).split('>')[0]
        #        if variable_type != 'CCVariable' and and variable_type != 'FCVariable' and variable_type != 'NCVariable':
        #            return None
        var_elem = ET.Element('DataArray')
        if number_type == 'long64':
          var_elem.attrib['type'] = 'UInt64'
        else:
          #var_elem.attrib['type'] = 'Float64'
          var_elem.attrib['type'] = 'Float32'
        # To add material number to variable: var_elem.attrib['Name'] = self.name + '_' + repr(self.index)
        #material = repr(self.index)
        var_elem.attrib['Name'] = self.name
        var_elem.attrib['NumberOfComponents'] = str(self.data.shape[1])
        var_elem.attrib['format'] = 'ascii'
        string_var_elem = str()
        for val in self.data.flatten():
          # set crazy small numbers to zero (grid variable number type double values)
          if abs(val)<tolerance:
            val = 0.0
          string_var_elem = string_var_elem + ' ' + repr(val)
          #print self.name +' ('+number_type+ ') mat'+material+'=' + repr(val) 
        var_elem.text = string_var_elem
        indent(var_elem)
        #print 'var_elem=',ET.dump(var_elem)
        return var_elem


# MAIN Subroutine 
if __name__ == "__main__":
  str_header=  '''
filename:   uda2vtk.py
  
summary:    Conversion utility to read Uintah *.uda output files and create VTK files
            for grid, particle and particle domain visualization and analysis

output:     VTK XML files are created in the <directory.uda> directory
            for reading with ParaView or other visualization tools

            grid.pvd            Uintah (traditional) grid files
            particles.pvd       Uintah (traditional) particle centroid files
            domains.pvd         Uintah particle domain files
  '''
  str_footer='''
comments:
The default is to write the particle domain VTK XML files.
More advanced visualization of particles domains requires additional variables 
to be saved in Uintah.  To visualize Convected Particle Domain Interpolator (CPDI)
or Convected Particle Tetrahedral/Triangular domain Interpolator (CPTI) particle domains,
the Uintah <p.scalefactor> variable needs to be saved during a simulation. 
  '''

  # timing
  time_start = time.time()

  # MPI
  if hasmpi:
    mpicomm = MPI.COMM_WORLD
  else:
    mpicomm = None

  # argument parsing
  parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=str_header,
    prog='python uda2vtk.py',
    epilog=str_footer)
  parser.add_argument(dest='udafile',metavar='<directory.uda>',action='store',default=True,
    help='the Uintah *.uda simulation output directory.')
  #parser.add_argument('-b','--binary',dest='do_binary',action='store_true',default=False,
  #  help='write binary VTK files.')
  parser.add_argument('-c','--compression',dest='do_compression',action='store_true',default=False,
    help='write compressed VTK files.')
  parser.add_argument('-d','--domains',dest='do_domains',action='store_false',default=True,
    help='do not write VTK particle domain files.')
  parser.add_argument('-g','--grid',dest='do_grid',action='store_true',default=False,
    help='write VTK grid files.')
  parser.add_argument('-j','--nprocs',dest='nprocs',action='store',default=1,
    help='run with multiple processors.')
  parser.add_argument('-p','--particles',dest='do_particles',action='store_true',default=False,
    help='write VTK particle centroid files.')
  parser.add_argument('-s','--size',dest='do_size',action='store_true',default=False,
    help='use p.size, p.deformationMeasure and grid spacing to write VTK particle domains.')
  parser.add_argument('-t','--timestep',dest='timedir',action='store',default=False,
    help='use a specific timestep directory for visualization.')
  parser.add_argument('-w','--wipe',dest='wipe_vtk',action='store_true',default=False,
    help='wipe VTK files from Uintah uda directory.')
  args = parser.parse_args()

  # read uintah uda directory and files
  if os.path.exists(args.udafile):  
    if args.wipe_vtk:
      os.system('rm -rf '+args.udafile+'/t*/*vtu')
      os.system('rm -rf '+args.udafile+'/t*/l*/*vtu')
      os.system('rm -rf '+args.udafile+'/t*/l*/*vti')
      os.system('rm -rf '+args.udafile+'/t*/*vti')
      os.system('rm -rf '+args.udafile+'/*pvd')
      print 'VTK files wiped from '+args.udafile
      sys.exit()
    else:
      UDAFILE=os.path.realpath(args.udafile)
      inUda=Uda(mpicomm=mpicomm)
      # perhaps find way to append additional timesteps or include timestep list?
      if args.timedir and os.path.exists(UDAFILE+'/'+args.timedir):
        inUda.read(UDAFILE,args.timedir)
      else:
        inUda.read(UDAFILE)
  else:
    parser.print_help()
    print '*** ERROR *** Uintah directory '+args.udafile+' is missing.'
    sys.exit()

  # write VTK output files
  print '\nWriting VTK output.'
  # grid files (grid.pvd, grid.vti)
  if args.nprocs:
    nprocs=int(args.nprocs)
  else:
    nprocs=1
  if args.do_grid:
    print '...grid output files started.'
    inUda.output_vtk_grid(UDAFILE)
    print '...grid output files complete.\n'
  # particle centroid files (particles.pvd, particles.pvtu, particles.m*.vtu)
  if args.do_particles:
    print '...particle centroid output files started.'
    inUda.output_vtk_particle(UDAFILE)
    print '...particle centroid output files complete.\n'
  # particle domain files (domains.pvd, domains.pvtu, domains.m*.vtu) 
  if args.do_domains:
    print '...particle domain output files started.'
    inUda.output_vtk_particle_domain(UDAFILE,nprocs)
    print '...particle domain output files complete.\n'
  # finished
  print '\nVTK output complete in %.4f s'%(time.time()-time_start)
