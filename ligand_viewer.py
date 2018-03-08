import os
import operator

# directory and files
# =============================================================================
# base directory for files
root_dir = '/directory_with_files'

# CSV file containing ligand code, correlation
csv_file = 'ligands.csv'

# functions for Coot
# =============================================================================
class ligand_viewer(object):

  def __init__(self, ligands=None, root_dir=None, mark_file='mark.csv',
               save_file='save.csv'):
    '''
    Parameters:
    ===========
    ligands   - list of tuples [ (ligand, correlation) ]
    root_dir  - base directory containing ligands
    mark_file - file listing marked ligands (output)
    save_file - file listing filenames of saved ligands (output)

    Notes:
    ======
    The default path for created files is the current working directory
    Marked ligands are unique
    Saved ligands are overwritten if they already exist
    Saved ligands will also be marked, but can be unmarked later
    Only 1 model exists, so imol (from Coot) is always 0

    Usage:
    ======
    Modify root_dir and csv_file so that this file can find the ligands
    Start Coot: coot --script ligand_viewer.py
    '[' in Coot will go back one ligand
    ']' in Coot will go forward one ligand
    '`' in Coot will mark the current ligand (tilde key)
    's' in Coot will save the current ligand coordinates
    '''

    # files
    self.root_dir = root_dir
    self.mark_file = mark_file
    self.save_file = save_file

    # data
    self.ligands = ligands
    self.marked_ligands = set()
    self.saved_ligands = set()

    self.current_ligand_index = 0

    # do not change map colors
    set_colour_map_rotation_for_map(0)

    # bind keys
    add_key_binding('Previous ligand', '[', self.prev_ligand)
    add_key_binding('Next ligand', ']', self.next_ligand)
    add_key_binding('Mark ligand', '`', self.mark_ligand)
    add_key_binding('Save ligand', 's', self.save_ligand)

    # prefill marked_ligands if mark_file already exists
    if (os.path.isfile(self.mark_file)):
      f = open(self.mark_file, 'r')
      lines = f.readlines()
      f.close()

      for line in lines:
        split_line = line.split(',')
        current_ligand = (split_line[0].strip(), split_line[1].strip())
        self.marked_ligands.add(current_ligand)

    # create new save_file
    f = open(self.save_file, 'w')
    f.close()

  # ---------------------------------------------------------------------------
  def coot_print(self, message):
    print '*' * 79
    print message
    print '*' * 79

  # ---------------------------------------------------------------------------
  def dir_exists(self, directory):
    if (os.path.isdir(directory)):
      return True
    else:
      self.coot_print('%s does not exist' % directory)
      return False

  # ---------------------------------------------------------------------------
  def file_exists(self, filename):
    if (os.path.isfile(filename)):
      return True
    else:
      self.coot_print('%s does not exist' % filename)
      return False

  # ---------------------------------------------------------------------------
  def show_ligand(self, ligand_index=0):

    # check index
    if ( (ligand_index < 0) or (ligand_index > (len(self.ligands) - 1)) ):
      print 'ligand_index out of range'
      if (ligand_index < 0):
        self.current_ligand_index = 0
      else:
        self.current_ligand_index = len(self.ligands) - 1
      return

    # check that ligand directory exists
    current_ligand = self.ligands[ligand_index]
    ligand = current_ligand[0]
    ligand_dir = os.path.join(self.root_dir,ligand)
    ligand_fit_dir = os.path.join(ligand_dir, 'LigandFit_run_1_')
    if (self.dir_exists((ligand_dir))):
      if (current_ligand in self.marked_ligands):
        message = '%i %s %s already marked'
      else:
        message = '%i %s %s'
      self.coot_print(message % (ligand_index, ligand, current_ligand[1]))
      model_file = os.path.join(ligand_fit_dir, 'ligand_fit_1.pdb')
      map_file = os.path.join(ligand_fit_dir, 'resolve_map.mtz')
      cif_file = os.path.join(ligand_fit_dir, ligand + '.sdf_ELBOW.cif')
      if (ligand.endswith('_sdf')):
        cif_file = os.path.join(ligand_fit_dir, 'elbow.' + ligand + '.001_ELBOW.cif')

      # read CIF
      if (self.file_exists(cif_file)):
        read_cif_dictionary(cif_file)

      # remove maps
      # for imol in map_molecule_list():
      #   close_molecule(imol)

      # display model and map
      if (self.file_exists(model_file)):
        if (len(molecule_number_list()) == 0):
          imol = handle_read_draw_molecule_with_recentre(model_file, 1)
          if (self.file_exists(map_file)):
            map_imol = auto_read_make_and_draw_maps(map_file)
        else:
          clear_and_update_model_molecule_from_file(0, model_file)

    # update location
    self.current_ligand_index = ligand_index

  # ---------------------------------------------------------------------------
  def next_ligand(self):
    self.show_ligand(self.current_ligand_index + 1)

  # ---------------------------------------------------------------------------
  def prev_ligand(self):
    self.show_ligand(self.current_ligand_index - 1)

  # ---------------------------------------------------------------------------
  def mark_ligand(self):
    current_ligand = self.ligands[self.current_ligand_index]
    if (current_ligand in self.marked_ligands):
      self.marked_ligands.remove(current_ligand)
      f = open(self.mark_file, 'w')
      for ligand in self.marked_ligands:
        f.write('%s, %s\n' % ligand)
      f.close()
      self.coot_print('%s unmarked' % current_ligand[0])
    else:
      self.marked_ligands.add(current_ligand)
      f = open(self.mark_file, 'a')
      f.write('%s, %s\n' % current_ligand)
      f.close()
      self.coot_print('%s marked' % current_ligand[0])

  # ---------------------------------------------------------------------------
  def save_ligand(self):
    current_ligand = self.ligands[self.current_ligand_index]
    if (current_ligand not in self.saved_ligands):
      self.saved_ligands.add(current_ligand)
      f = open(self.save_file, 'a')
      f.write('%s, %s\n' % current_ligand)
      f.close()
    if (current_ligand not in self.marked_ligands):
      self.marked_ligands.add(current_ligand)
      f = open(self.mark_file, 'a')
      f.write('%s, %s\n' % current_ligand)
      f.close()

    filename = current_ligand[0] + '_new.pdb'
    save_coordinates(0, filename)

    self.coot_print('%s saved' % current_ligand[0])

# initial setup
# =============================================================================

f = open(csv_file, 'r')
lines = f.readlines()
f.close()

# parse ligand data
ligands = list()
for line in lines:
  split_line = line.split(',')
  current_ligand = split_line[0].strip()
  correlation = split_line[1].strip()

  # isolate M2MDBxxxxxx
  split_line = current_ligand.split('.')
  if ('M2MDB' in split_line[0]):
    current_ligand = split_line[0]
  else:
    current_ligand = split_line[1]

  ligands.append((current_ligand, correlation))

# sort by correlation (descending)
ligands.sort(key=operator.itemgetter(1), reverse=True)

# initialize class
lv = ligand_viewer(ligands=ligands,root_dir=root_dir)

# load first ligand
lv.show_ligand()

# =============================================================================
# end of file
