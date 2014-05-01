#=============================================================================================
# Set up hydrated and vacuum simulations of small molecule in water.
#
# PROTOCOL
#
# * Construct molecule from IUPAC name (protonation and tautomeric states are heuristically guessed) [OpenEye OEMol]
# * Generate multiple likely conformations to use for replicates [OpenEye Omega]
# * Parameterize the molecule with GAFF [AmberTools Antechamber]
# * Set up the solvated system [AmberTools LEaP]
# 
# Written by John D. Chodera <jchodera@gmail.com>
#
# TODO:
# - Write out sander mdin files for minimization and dynamics.
#=============================================================================================

#=============================================================================================
# Imports
#=============================================================================================

import commands
import os    
from numpy import *

from openeye.oechem import *
from openeye.oeomega import *
from openeye.oeiupac import *
from openeye.oeshape import *
try:
   from openeye.oequacpac import * #DLM added 2/25/09 for OETyperMolFunction; replacing oeproton
except:
   from openeye.oeproton import * #GJR temporary fix because of old version of openeye tools
from openeye.oeiupac import *
from openeye.oeszybki import *

#=============================================================================================
# PARAMETERS
#=============================================================================================

clearance = 9.0 # clearance around solute for box construction, in Angstroms

sander_minimization_mdin = """\
initial minimization prior to MD
 &cntrl
  imin   = 1,
  maxcyc = 500,
  ncyc   = 250,
  ntb    = 0,
  igb    = 0,
  cut    = 9
 /
"""

sander_dynamics_mdin = """\
initial minimization prior to MD
 &cntrl
  imin   = 0,
  maxcyc = 5000,
  ntb    = 0,
  igb    = 0,
  cut    = 9
 /
"""

#=============================================================================================
# SUBROUTINES
#=============================================================================================

def createMoleculeFromIUPAC(name, verbose = False, charge = None, strictTyping = None, strictStereo = True):
   """Generate a small molecule from its IUPAC name.

   ARGUMENTS
     IUPAC_name (string) - IUPAC name of molecule to generate

   OPTIONAL ARGUMENTS
     verbose (boolean) - if True, subprocess output is shown (default: False)
     charge (int) - if specified, a form of this molecule with the desired charge state will be produced (default: None)
     strictTyping (boolean) -- if set, passes specified value to omega (see documentation for expandConformations)
     strictStereo (boolean) -- if True, require stereochemistry to be specified before running omega. If False, don't (pick random stereoisomer if not specified). If not specified (None), do whatever omega does by default (varies with version). Default: True.

   RETURNS
     molecule (OEMol) - the molecule

   NOTES
     OpenEye LexiChem's OEParseIUPACName is used to generate the molecle.
     The molecule is normalized by adding hydrogens.
     Omega is used to generate a single conformation.
     Also note that atom names will be blank coming from this molecule. They are assigned when the molecule is written, or one can assign using OETriposAtomNames for example.

   EXAMPLES
     # Generate a mol2 file for phenol.
     molecule = createMoleculeFromIUPAC('phenol')
     
   """

   # Create an OEMol molecule from IUPAC name.
   molecule = OEMol() # create a molecule
   status = OEParseIUPACName(molecule, name) # populate the molecule from the IUPAC name

   # Normalize the molecule.
   normalizeMolecule(molecule)

   # Generate a conformation with Omega
   omega = OEOmega()
   if strictStereo<>None:
        omega.SetStrictStereo(strictStereo)

   #omega.SetVerbose(verbose)
   #DLM 2/27/09: Seems to be obsolete in current OEOmega
   if strictTyping != None:
     omega.SetStrictAtomTypes( strictTyping) 
   
   omega.SetIncludeInput(False) # don't include input
   omega.SetMaxConfs(1) # set maximum number of conformations to 1
   omega(molecule) # generate conformation      

   if (charge != None):
      # Enumerate protonation states and select desired state.
      protonation_states = enumerateStates(molecule, enumerate = "protonation", verbose = verbose)
      for molecule in protonation_states:
         if formalCharge(molecule) == charge:
            # Return the molecule if we've found one in the desired protonation state.
            return molecule
      if formalCharge(molecule) != charge:
         print "enumerateStates did not enumerate a molecule with desired formal charge."
         print "Options are:"
         for molecule in protonation_states:
            print "%s, formal charge %d" % (molecule.GetTitle(), formalCharge(molecule))
         raise "Could not find desired formal charge."
    
   # Return the molecule.
   return molecule

def writeMolecule(molecule, filename, substructure_name = 'MOL', preserve_atomtypes = False):
   """Write a molecule to a file in any format OpenEye autodetects from filename (such as .mol2).
   WARNING: The molecule will be standardized before writing by the high-level OEWriteMolecule function.
   OEWriteConstMolecule is used, to avoid changing the molecule you pass in.

   ARGUMENTS
     molecule (OEMol) - the molecule to be written
     filename (string) - the file to write the molecule to (type autodetected from filename)

   OPTIONAL ARGUMENTS
     substructure_name (String) - if a mol2 file is written, this is used for the substructure name (default: 'MOL')
     preserve_atomtypes (bool) - if True, a mol2 file will be written with atom types preserved

   RETURNS
     None

   NOTES
     Multiple conformers are written.

   EXAMPLES
     # create a molecule
     molecule = createMoleculeFromIUPAC('phenol')
     # write it as a mol2 file
     writeMolecule(molecule, 'phenol.mol2')
   """

   # Open output stream.
   ostream = oemolostream()
   ostream.open(filename)

   # Define internal function for writing multiple conformers to an output stream.
   def write_all_conformers(ostream, molecule):
      # write all conformers of each molecule
      for conformer in molecule.GetConfs():
         if preserve_atomtypes: OEWriteMol2File(ostream, conformer)
         else: OEWriteConstMolecule(ostream, conformer)
      return

   # If 'molecule' is actually a list of molecules, write them all.
   if type(molecule) == type(list()):
      for individual_molecule in molecule:
         write_all_conformers(ostream, individual_molecule)
   else:
      write_all_conformers(ostream, molecule)

   # Close the stream.
   ostream.close()

   # Replace substructure name if mol2 file.
   suffix = os.path.splitext(filename)[-1]
   if (suffix == '.mol2' and substructure_name != None):
      modifySubstructureName(filename, substructure_name)

   return

def readMolecule(filename, normalize = False):
   """Read in a molecule from a file (such as .mol2).

   ARGUMENTS
     filename (string) - the name of the file containing a molecule, in a format that OpenEye autodetects (such as .mol2)

   OPTIONAL ARGUMENTS
     normalize (boolean) - if True, molecule is normalized (renamed, aromaticity, protonated) after reading (default: False)

   RETURNS
     molecule (OEMol) - OEMol representation of molecule

   EXAMPLES
     # read a mol2 file
     molecule = readMolecule('phenol.mol2')
     # works with any type of file that OpenEye autodetects
     molecule = readMolecule('phenol.sdf')
   """

   # Open input stream.
   istream = oemolistream()
   istream.open(filename)

   # Create molecule.
   molecule = OEMol()   

   # Read the molecule.
   OEReadMolecule(istream, molecule)

   # Close the stream.
   istream.close()

   # Normalize if desired.
   if normalize: normalizeMolecule(molecule)

   return molecule


def enumerateStates(molecules, enumerate = "protonation", consider_aromaticity = True, maxstates = 200, verbose = True):
    """Enumerate protonation or tautomer states for a list of molecules.

    ARGUMENTS
      molecules (OEMol or list of OEMol) - molecules for which states are to be enumerated

    OPTIONAL ARGUMENTS
      enumerate - type of states to expand -- 'protonation' or 'tautomer' (default: 'protonation')
      verbose - if True, will print out debug output

    RETURNS
      states (list of OEMol) - molecules in different protonation or tautomeric states

    TODO
      Modify to use a single molecule or a list of molecules as input.
      Apply some regularization to molecule before enumerating states?
      Pick the most likely state?
      Add more optional arguments to control behavior.
    """

    # If 'molecules' is not a list, promote it to a list.
    if type(molecules) != type(list()):
       molecules = [molecules]

    # Check input arguments.
    if not ((enumerate == "protonation") or (enumerate == "tautomer")):
        raise "'enumerate' argument must be either 'protonation' or 'tautomer' -- instead got '%s'" % enumerate

    # Create an internal output stream to expand states into.
    ostream = oemolostream()
    ostream.openstring()
    ostream.SetFormat(OEFormat_SDF)
    
    # Default parameters.
    only_count_states = False # enumerate states, don't just count them

    # Enumerate states for each molecule in the input list.
    states_enumerated = 0
    for molecule in molecules:
        if (verbose): print "Enumerating states for molecule %s." % molecule.GetTitle()
        
        # Dump enumerated states to output stream (ostream).
        if (enumerate == "protonation"): 
            # Create a functor associated with the output stream.
            functor = OETyperMolFunction(ostream, consider_aromaticity, False, maxstates)
            # Enumerate protonation states.
            if (verbose): print "Enumerating protonation states..."
            states_enumerated += OEEnumerateFormalCharges(molecule, functor, verbose)        
        elif (enumerate == "tautomer"):
            # Create a functor associated with the output stream.
            functor = OETautomerMolFunction(ostream, consider_aromaticity, False, maxstates)
            # Enumerate tautomeric states.
            if (verbose): print "Enumerating tautomer states..."
            states_enumerated += OEEnumerateTautomers(molecule, functor, verbose)    
    print "Enumerated a total of %d states." % states_enumerated

    # Collect molecules from output stream into a list.
    states = list()
    if (states_enumerated > 0):    
        state = OEMol()
        istream = oemolistream()
        istream.openstring(ostream.GetString())
        istream.SetFormat(OEFormat_SDF)
        while OEReadMolecule(istream, state):
           states.append(OEMol(state)) # append a copy

    # Return the list of expanded states as a Python list of OEMol() molecules.
    return states

def formalCharge(molecule):
   """Report the net formal charge of a molecule.

   ARGUMENTS
     molecule (OEMol) - the molecule whose formal charge is to be determined

   RETURN VALUES
     formal_charge (integer) - the net formal charge of the molecule

   EXAMPLE
     net_charge = formalCharge(molecule)
   """

   # Create a copy of the molecule.
   molecule_copy = OEMol(molecule)

   # Assign formal charges.
   OEFormalPartialCharges(molecule_copy)

   # Compute net formal charge.
   formal_charge = int(round(OENetCharge(molecule_copy)))

   # return formal charge
   return formal_charge

def normalizeMolecule(molecule):
   """Normalize the molecule by checking aromaticity, adding explicit hydrogens, and renaming by IUPAC name.

   ARGUMENTS
     molecule (OEMol) - the molecule to be normalized.

   EXAMPLES
     # read a partial molecule and normalize it
     molecule = readMolecule('molecule.sdf')
     normalizeMolecule(molecule)
   """
   
   # Find ring atoms and bonds
   # OEFindRingAtomsAndBonds(molecule) 
   
   # Assign aromaticity.
   OEAssignAromaticFlags(molecule, OEAroModelOpenEye)   

   # Add hydrogens.
   OEAddExplicitHydrogens(molecule)

   # Set title to IUPAC name.
   name = OECreateIUPACName(molecule)
   molecule.SetTitle(name)

   return molecule

def expandConformations(molecule, maxconfs = None, threshold = None, include_original = False, torsionlib = None, verbose = False, strictTyping = None, strictStereo = None):   
   """Enumerate conformations of the molecule with OpenEye's Omega after normalizing molecule. 

   ARGUMENTS
   molecule (OEMol) - molecule to enumerate conformations for

   OPTIONAL ARGUMENTS
     include_original (boolean) - if True, original conformation is included (default: False)
     maxconfs (integer) - if set to an integer, limits the maximum number of conformations to generated -- maximum of 120 (default: None)
     threshold (real) - threshold in RMSD (in Angstroms) for retaining conformers -- lower thresholds retain more conformers (default: None)
     torsionlib (string) - if a path to an Omega torsion library is given, this will be used instead (default: None)
     verbose (boolean) - if True, omega will print extra information
     strictTyping (boolean) -- if specified, pass option to SetStrictAtomTypes for Omega to control whether related MMFF types are allowed to be substituted for exact matches.
     strictStereo (boolean) -- if specified, pass option to SetStrictStereo; otherwise use default.

   RETURN VALUES
     expanded_molecule - molecule with expanded conformations

   EXAMPLES
     # create a new molecule with Omega-expanded conformations
     expanded_molecule = expandConformations(molecule)

     
   """
   # Initialize omega
   omega = OEOmega()
   if strictTyping != None:
     omega.SetStrictAtomTypes( strictTyping)
   if strictStereo != None:
     omega.SetStrictStereo( strictStereo )
   #Set atom typing options

   # Set verbosity.
   #omega.SetVerbose(verbose)
   #DLM 2/27/09: Seems to be obsolete in current OEOmega

   # Set maximum number of conformers.
   if maxconfs:
      omega.SetMaxConfs(maxconfs)
     
   # Set whether given conformer is to be included.
   omega.SetIncludeInput(include_original)
   
   # Set RMSD threshold for retaining conformations.
   if threshold:
      omega.SetRMSThreshold(threshold) 
 
   # If desired, do a torsion drive.
   if torsionlib:
      omega.SetTorsionLibrary(torsionlib)

   # Create copy of molecule.
   expanded_molecule = OEMol(molecule)   

   # Enumerate conformations.
   omega(expanded_molecule)


   # verbose output
   if verbose: print "%d conformation(s) produced." % expanded_molecule.NumConfs()

   # return conformationally-expanded molecule
   return expanded_molecule

def modifySubstructureName(mol2file, name):
   """Replace the substructure name (subst_name) in a mol2 file.

   ARGUMENTS
     mol2file (string) - name of the mol2 file to modify
     name (string) - new substructure name

   NOTES
     This is useful becuase the OpenEye tools leave this name set to <0>.
     The transformation is only applied to the first molecule in the mol2 file.

   TODO
     This function is still difficult to read.  It should be rewritten to be comprehensible by humans.
   """

   # Read mol2 file.
   file = open(mol2file, 'r')
   text = file.readlines()
   file.close()

   # Find the atom records.
   atomsec = []
   ct = 0
   while text[ct].find('<TRIPOS>ATOM')==-1:
     ct+=1
   ct+=1
   atomstart = ct
   while text[ct].find('<TRIPOS>BOND')==-1:
     ct+=1
   atomend = ct

   atomsec = text[atomstart:atomend]
   outtext=text[0:atomstart]
   repltext = atomsec[0].split()[7] # mol2 file uses space delimited, not fixed-width

   # Replace substructure name.
   for line in atomsec:
     # If we blindly search and replace, we'll tend to clobber stuff, as the subst_name might be "1" or something lame like that that will occur all over. 
     # If it only occurs once, just replace it.
     if line.count(repltext)==1:
       outtext.append( line.replace(repltext, name) )
     else:
       # Otherwise grab the string left and right of the subst_name and sandwich the new subst_name in between. This can probably be done easier in Python 2.5 with partition, but 2.4 is still used someplaces.
       # Loop through the line and tag locations of every non-space entry
       blockstart=[]
       ct=0
       c=' '
       for ct in range(len(line)):
         lastc = c
         c = line[ct]
         if lastc.isspace() and not c.isspace():
           blockstart.append(ct)
       line = line[0:blockstart[7]] + line[blockstart[7]:].replace(repltext, name, 1)
       outtext.append(line)
       
   # Append rest of file.
   for line in text[atomend:]:
     outtext.append(line)
     
   # Write out modified mol2 file, overwriting old one.
   file = open(mol2file,'w')
   file.writelines(outtext)
   file.close()

   return

def write_file(filename, contents):
    """Write the specified contents to a file.

    ARGUMENTS
      filename (string) - the file to be written
      contents (string) - the contents of the file to be written

    """

    outfile = open(filename, 'w')

    if type(contents) == list:
        for line in contents:
            outfile.write(line)
    elif type(contents) == str:
        outfile.write(contents)
    else:
        raise "Type for 'contents' not supported: " + repr(type(contents))

    outfile.close()

    return

def read_file(filename):
    """Read contents of the specified file.

    ARGUMENTS
      filename (string) - the name of the file to be read

    RETURNS
      lines (list of strings) - the contents of the file, split by line
    """

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    return lines
#=============================================================================================
def setupHydrationCalculation(solute, verbose=True):
    """Set up an absolute alchemical hydration free energy calculation for the given molecule.

    ARGUMENTS
      solute (OEMol) - the molecule for which hydration free energy is to be computed (with fully explicit hydrogens) in the desired protonation state.

    OPTIONAL ARGUMENTS
      verbose (boolean) - if True, extra debug information will be printed (default: True)

    NOTES
      A directory will be created 'molecules/[molecule name]' as obtained from molecule.GetTitle().
    
    """

    # get current directory
    current_path = os.getcwd()

    # Center the solute molecule.
    OECenter(solute)

    # get molecule name
    solute_name = molecule.GetTitle()
    if verbose: print solute_name

    # create molecule path/directory
    work_path = os.path.abspath(os.path.join('molecules', solute_name))
    os.makedirs(work_path)

    # SET UP SOLUTE TOPOLOGY

    if verbose: print "\nCONSTRUCTING SOLUTE TOPOLOGY"

    # Get formal charge of ligand.
    solute_charge = formalCharge(solute)
    if verbose: print "solute formal charge is %d" % solute_charge
    
    # Write molecule with explicit hydrogens to mol2 file.
    print "Writing solute mol2 file..."
    solute_mol2_filename = os.path.abspath(os.path.join(work_path, 'solute.mol2'))
    writeMolecule(solute, solute_mol2_filename)

    # Set substructure name (which will become residue name).
    print "Modifying molecule name..."
    modifySubstructureName(solute_mol2_filename, 'MOL')

    # Run antechamber to assign GAFF atom types.
    print "Running antechamber..."
    os.chdir(work_path)
    gaff_mol2_filename = os.path.join(work_path, 'solute.gaff.mol2')
    charge_model = 'bcc'
    command = 'antechamber -i %(solute_mol2_filename)s -fi mol2 -o solute.gaff.mol2 -fo mol2 -c %(charge_model)s -nc %(solute_charge)d > antechamber.out' % vars()
    if verbose: print command
    output = commands.getoutput(command)
    if verbose: print output
    os.chdir(current_path)

    # Generate frcmod file for additional GAFF parameters.
    solute_frcmod_filename = os.path.join(work_path, 'frcmod.solute')
    command = 'parmchk -i %(gaff_mol2_filename)s -f mol2 -o %(solute_frcmod_filename)s' % vars()
    if verbose: print command
    output = commands.getoutput(command)
    if verbose: print output

    # Run LEaP to generate topology / coordinates.
    solute_prmtop_filename = os.path.join(work_path,'solute.prmtop')
    solute_crd_filename = os.path.join(work_path,'solute.inpcrd')
    solute_off_filename = os.path.join(work_path, 'solute.off')
    
    tleap_input_filename = os.path.join(work_path, 'setup-solute.leap.in')
    tleap_output_filename = os.path.join(work_path, 'setup-solute.leap.out')
    contents = """
# Load GAFF parameters.
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(solute_frcmod_filename)s

# load solute
solute = loadMol2 %(gaff_mol2_filename)s

# check the solute
check solute

# report net charge
charge solute

# save AMBER parameters
saveAmberParm solute %(solute_prmtop_filename)s %(solute_crd_filename)s

# write .off file
saveOff solute %(solute_off_filename)s

# exit
quit
""" % vars()
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)

    # extract total charge
    solute_charge = commands.getoutput('grep "Total unperturbed charge" %(tleap_output_filename)s | cut -c 27-' % vars())
    solute_charge = int(round(float(solute_charge))) # round to nearest whole charge
    if verbose: print "solute charge is %d" % solute_charge    

    # PREPARE SOLVATED SOLUTE

    print "\nPREPARING SOLVATED SOLUTE"

    # create the directory if it doesn't exist
    solvent_path = os.path.join(work_path, 'solvent')
    if not os.path.exists(solvent_path):
        os.makedirs(solvent_path)

    # solvate the solute
    print "Solvating the solute with tleap..."
    system_prmtop_filename = os.path.join(solvent_path,'system.prmtop')
    system_crd_filename = os.path.join(solvent_path,'system.inpcrd')
    tleap_input_filename = os.path.join(solvent_path, 'setup-system.leap.in')
    tleap_output_filename = os.path.join(solvent_path, 'setup-system.leap.out')
    clearance = globals()['clearance'] # clearance around solute (in A)
    contents = """
source leaprc.ff10
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(solute_frcmod_filename)s

# Load solute.
loadOff %(solute_off_filename)s

# Create system.
system = combine { solute }
""" % vars()
    # add counterions
    if (solute_charge != 0):
        nions = abs(solute_charge)
        if solute_charge < 0: iontype = 'Na+'
        if solute_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions.
addions system %(iontype)s %(nions)d
""" % vars()
    #
    contents += """
# Solvate in truncated octahedral box.
solvateBox system TIP3PBOX %(clearance)f iso

# Check the system
check system

# Write the system
saveamberparm system %(system_prmtop_filename)s %(system_crd_filename)s
""" % vars()    
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)    

    # make a PDB file for checking
    print "Converting system to PDB..."
    os.chdir(solvent_path)
    command = 'cat system.inpcrd | ambpdb -p system.prmtop > system.pdb' % vars()
    output = commands.getoutput(command)
    print output
    os.chdir(current_path)

    # SET UP VACUUM SIMULATION

    # construct pathname for vacuum simulations
    vacuum_path = os.path.join(work_path, 'vacuum')
    if not os.path.exists(vacuum_path):
        os.makedirs(vacuum_path)

    # solvate the solute
    print "Preparing vacuum solute with tleap..."
    system_prmtop_filename = os.path.join(vacuum_path,'system.prmtop')
    system_crd_filename = os.path.join(vacuum_path,'system.inpcrd')
    tleap_input_filename = os.path.join(vacuum_path, 'setup-system.leap.in')
    tleap_output_filename = os.path.join(vacuum_path, 'setup-system.leap.out')
    clearance = 50.0 # clearance in A
    contents = """
source leaprc.ff10
source leaprc.gaff

# load antechamber-generated additional parameters
mods = loadAmberParams %(solute_frcmod_filename)s

# Load solute.
loadOff %(solute_off_filename)s

# Create system.
system = combine { solute }
""" % vars()
    # add counterions
    if (solute_charge != 0):
        nions = abs(solute_charge)
        if solute_charge < 0: iontype = 'Na+'
        if solute_charge > 0: iontype = 'Cl-'
        contents += """
# Add counterions.
addions system %(iontype)s %(nions)d
""" % vars()
    #
    contents += """

# Create big box.
setBox system centers %(clearance)f iso

# Check the system
check system

# Write the system
saveamberparm system %(system_prmtop_filename)s %(system_crd_filename)s
""" % vars()    
    write_file(tleap_input_filename, contents)
    command = 'tleap -f %(tleap_input_filename)s > %(tleap_output_filename)s' % vars()
    output = commands.getoutput(command)    
    
    # make a PDB file for checking
    print "Converting system to PDB..."
    os.chdir(vacuum_path)
    command = 'cat system.inpcrd | ambpdb -p system.prmtop > system.pdb' % vars()
    output = commands.getoutput(command)
    print output
    os.chdir(current_path)

    return


#=============================================================================================
# MAIN
#=============================================================================================

# List of tuples describing molecules to construct:
# ('common name to use for directory', 'IUPAC name', formal_charge to select)
molecules = [
    ('methane', 'methane', 0), # methane
    ('benzene', 'benzene', 0), # benzene
    ('ethanol', 'ethanol', 0), # ethanol
    ('phenol', 'phenol', 0), # phenol
    ('N-methylacetamide', 'N-methylacetamide', 0),
]

for (molecule_name, molecule_iupac_name, formal_charge) in molecules:
    print "Setting up %s" % molecule_name

    # Create the molecule from the common name with the desired formal charge.
    molecule = createMoleculeFromIUPAC(molecule_iupac_name, charge = formal_charge)
    if not molecule:
        print "Could not build '%(molecule_name)s' from IUPAC name '%(molecule_iupac_name)s', skipping..." % vars()
        continue
    print "Net charge is %d" % formalCharge(molecule)

    # Replace the title with the common name
    molecule.SetTitle(molecule_name)

    # Expand set of conformations so we have multiple conformations to start from.
    expandConformations(molecule)
    
    # Set up multiple replicates of the hydration free energy calculation.
    setupHydrationCalculation(molecule)

    print "\n\n"

