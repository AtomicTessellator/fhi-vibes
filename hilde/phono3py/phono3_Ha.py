from .phono import PhononsBase
import numpy              as     np
from copy import deepcopy as copy

class Phonons3(PhononsBase):
    def __init__(self, setup, structure):
        import os
        from   HIGHaims.structure import Cell
        from   HIGHaims.tasks     import Forces
        super(Phonons3, self).__init__(setup, structure)
        # Phonopy imports
        from phono3py.phonon3 import Phono3py

        if setup.fc3.symmetrize_fc:
            self.symmetrize_fc = True
        else:
            self.symmetrize_fc = False

        phonon = Phono3py(self.phonopy_unitcell,
                          supercell_matrix        = self.fc3_sc_matrix,
                          phonon_supercell_matrix = self.fc2_sc_matrix,
                          mesh                    = self.qmesh,
                          symprec                 = setup.structure.symprec,
                          is_symmetry             = True,
                          log_level               = 2,
                          frequency_factor_to_THz = self.const.eVtoTHz)

        # Primitive and supercell
        fc3_supercell        = Cell(phonon.get_supercell())
        fc2_supercell        = Cell(phonon.get_phonon_supercell())
        primitive            = Cell(self.phonopy_unitcell)

        # generate displacements
        if setup.fc3.cutoff:
          self.fc3_cutoff  = setup.fc3.cutoff
        else:
          self.fc3_cutoff  = None
        phonon.generate_displacements(distance             = setup.fc3.finite_disp,
                                      cutoff_pair_distance = self.fc3_cutoff,
                                      is_plusminus         ='auto',
                                      is_diagonal          = True)
        #
        self.fc3_cells_with_disps = phonon.get_supercells_with_displacements()
        self.fc2_cells_with_disps = phonon.get_phonon_supercells_with_displacements()

        fc3_tasks         = []
        fc3_workflow_dirs_rel = []
        for ii, cell in enumerate(self.fc3_cells_with_disps):
            if not cell:
                continue
            struc   = Cell(cell, symprec = None)
            wd_temp = os.path.join('fc3_displacements', 'fc3_displaced_cell_{:08d}'.format(ii))
            fc3_tasks.append(Forces({**setup.dft, **setup.fc3}, struc))
            fc3_workflow_dirs_rel.append(wd_temp)
        #
        # fc2 tasks and subdirectories
        fc2_tasks         = []
        fc2_workflow_dirs_rel = []
        #
        for ii, cell in enumerate(self.fc2_cells_with_disps):
            struc   = Cell(cell, symprec = None)
            wd_temp = os.path.join('fc2_displacements', 'fc2_displaced_cell_{:08d}'.format(ii))
            fc2_tasks.append(Forces({**setup.dft, **setup.fc2}, struc))
            fc2_workflow_dirs_rel.append(wd_temp)
        #
        self.fc2_tasks         = fc2_tasks

        # fc3 workflow parent directory
        workflow_dir = os.path.join(
            f'{setup.database.workdir_workflow}',
            f'{self.get_name()}'
        )

        # count the number of tasks:
        self.ntasks = [len(fc2_tasks), len(fc3_tasks)]

        # Write setup objects
        self.fc3_handler.setup.write(os.path.join(workflow_dir, 'setup.log'))

        # Write disp_fc3.yaml:
        from phono3py.file_IO import write_disp_fc3_yaml
        dds = phonon.get_displacement_dataset()
        num_disps, num_disp_files = write_disp_fc3_yaml(dds, fc3_supercell,
            filename = os.path.join(workflow_dir, 'fc3_disp.yaml'))

        # attach the phonon object to the workflow
        self.phonon         = phonon
        self.workflow_dir   = workflow_dir
        self.fc2_setup      = fc2_setup

#
    # Naming
    def get_hash(self, short = True):
        """ Hash is build from the hashes of the
        primitive cell,
        supercell,
        control for the tasks,
        first supercell with displacement (to be sure)
        """
        pc          = self.primitive
        sc2         = self.fc2_supercell
        sc3         = self.fc3_supercell
        cutoff      = len(self.fc3_tasks)
        task2       = self.fc2_tasks[0]
        task3       = self.fc3_tasks[0]
        phono3_hash = '_'.join(
            [pc.get_hash(short),
             sc2.get_hash(short),
             sc3.get_hash(short),
             task2.control.get_hash(short),
             f'{cutoff}',
             task3.control.get_hash(short),
             task3.structure.get_hash(short)]
            )
        return phono3_hash
    #
    def get_name(self, short = True):
        """ Some more human readable extension to the bare hash"""
        pc          = self.primitive
        xc          = self.setup.xc
        bset        = self.setup.basisset
        smf2        = self.fc2_sc_matrix.flatten()
        smf3        = self.fc3_sc_matrix.flatten()
        cutoff      = len(self.fc3_tasks)
        dist        = self.setup.phonopy_dist
        vol         = pc.get_volume()
        workflow    = 'phono3'
        phonon_hash = self.get_hash(short)
        workflow_name = '_'.join(
            [workflow,
             pc.sysname,
             f'{xc}',
             f'{bset}',
             f'{vol:.2f}',
             f'{cutoff}',
             f'{dist:.3f}',
             '_'.join(
                 [f'{smf2[ii]}{smf2[ii+1]}{smf2[ii+2]}' for ii in range(0, 9, 3)]),
             '_'.join(
                 [f'{smf3[ii]}{smf3[ii+1]}{smf3[ii+2]}' for ii in range(0, 9, 3)]),
            phonon_hash]
            )
        return workflow_name

    # Execute fc2
    def execute_fc2(self):
        # produce and store if not already calculated:
        import os
        from phono3py.file_IO import (read_fc2_from_hdf5, write_fc2_to_hdf5)
        phonon          = self.phonon
        if self.symmetrize_fc:
            fc2_fname = os.path.join(self.workflow_dir, 'fc2.symmetrized.hdf5')
        else:
            fc2_fname = os.path.join(self.workflow_dir, 'fc2.non_symmetrized.hdf5')
        if (not os.path.exists(fc2_fname)) :
            # If force constants do not exists yet, make force calculations
            fc2_status       = self.fc2_handler.execute()
            fc2_forces       = self.fc2_handler.get_result()

            # produce force constants
            if self.symmetrize_fc:
                phonon.produce_fc2(
                        fc2_forces,
                        is_translational_symmetry   = True,
                        is_permutation_symmetry     = True
                        )
            else:
                phonon.produce_fc2(fc2_forces)

            fc2 = phonon.get_fc2()
            write_fc2_to_hdf5(fc2, filename = fc2_fname)
        else:
            fc2 = read_fc2_from_hdf5(filename = fc2_fname)
            phonon.set_fc2(fc2)

        return 1
    #
    # execute fc2 and fc3
    def execute_fc3(self):
        # produce and store if not already calculated:
        import os
        from phono3py.file_IO import (read_fc3_from_hdf5, write_fc3_to_hdf5)
        phonon          = self.phonon
        if self.symmetrize_fc:
            fc3_fname = os.path.join(self.workflow_dir, 'fc3.symmetrized.hdf5')
        else:
            fc3_fname = os.path.join(self.workflow_dir, 'fc3.non_symmetrized.hdf5')
        if (not os.path.exists(fc3_fname)) :
            # if force constants do not exists, calculate forces and produce them
            fc3_status       = self.fc3_handler.execute()
            fc3_forces       = self.fc3_handler.get_result()
            # restore shape of fc3_forces if cutoff is used
            forces_new = []
            zero_force = np.zeros(fc3_forces[0].shape)
            counter    = 0
            for ii, cell in enumerate(self.fc3_cells_with_disps):
                if cell == None:
                    forces_new.append(zero_force)
                else:
                    forces_new.append(fc3_forces[counter])
                    counter += 1
            #
            fc3_forces = forces_new

            if self.symmetrize_fc:
                phonon.produce_fc3(
                        fc3_forces,
                        is_translational_symmetry   = True,
                        is_permutation_symmetry     = True,
                        is_permutation_symmetry_fc2 = True
                        )
            else:
                phonon.produce_fc3(fc3_forces)
            #
            fc3 = phonon.get_fc3()
            write_fc3_to_hdf5(fc3, filename = fc3_fname)



        else:
            fc3 = read_fc3_from_hdf5(filename = fc3_fname)
            phonon.set_fc3(fc3)
        #
        return 1
    # execute all
    def execute(self):
        ex1 = self.execute_fc2()
        ex2 = self.execute_fc3()
        if ex1 == 1 and ex2 == 1:
            return 1
        else:
            return -1

    def get_result(self):
        return self.phonon

    # Greeter
    def greet(self):
        print(f'\n[Ha] phono3py workflow. These are the settings I see:')
        print(f'  Material:                {self.primitive.sysname}')
        print(f'  XC functional:           {self.setup.xc}\n')
        print(f'  sc_acc_rho:              {self.setup.acc.rho}')
        print(f'  sc_acc_forces:           {self.setup.acc.forces}\n')
        print(f'  FC2 Supercell matrix:    {self.fc2_sc_matrix.flatten()}')
        print(f'  (number of atoms in sc): {self.fc2_supercell.Natoms}')
        print(f'  (number of disp.s)       {len(self.fc2_cells_with_disps)}')
        print(f'  (kpoints):               {self.fc2_setup.kgrid}\n')
        print(f'  FC3 Supercell matrix:    {self.fc3_sc_matrix.flatten()}')
        print(f'  FC3 cutoff:              {self.fc3_cutoff}')
        print(f'  FC2 displacement:        {self.setup.phonopy_dist}')
        print(f'  FC3 displacement:        {self.setup.phono3py_dist}')
        print(f'  (number of atoms in sc): {self.fc3_supercell.Natoms}')
        print(f'  (# disps)                {len([cc for cc in self.fc3_cells_with_disps if cc])}')
        print(f'  (# disps w/o cutoff)     {len(self.fc3_cells_with_disps)}')
        print(f'  (kpoints):               {self.setup.kgrid}\n')
        print(f'  Symmetrize FC:           {self.symmetrize_fc}\n')
        print(f'  Q-mesh for post process: {self.setup.qmesh}')
        print(f'  Include isotope effect:  {self.setup.kappa_isotope}\n')
        print(f'Identifier: {self.get_name()}')
        print(f'  -> {self.workflow_dir}')
