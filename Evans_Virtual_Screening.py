{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "a48cbb6d-41f3-4ca5-aed9-95eea8032b83",
   "metadata": {},
   "outputs": [],
   "source": [
    "import os\n",
    "os.chdir(path = r\"C:\\Users\\Admin\\Desktop\\VIRTUAL SCREENING\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "b72c824c-eaad-4419-b41e-3818e3915edd",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "'C:\\\\Users\\\\Admin\\\\Desktop\\\\VIRTUAL SCREENING'"
      ]
     },
     "execution_count": 2,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "os.getcwd()\n"
   ]
  },
  {
   "cell_type": "raw",
   "id": "d0f4a335-ce22-4bc3-ab4e-dd9d183acb0c",
   "metadata": {},
   "source": [
    "Ligand Preparation Step"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "58c15208-44ee-4402-a47a-86f493541bbe",
   "metadata": {},
   "outputs": [],
   "source": [
    "import os\n",
    "import subprocess\n",
    "\n",
    "# Define the input file name and output Folder\n",
    "input_file = 'query_results_zinc_2.sdf'\n",
    "output_dir = 'Processed_Ligand'\n",
    "os.makedirs(output_dir, exist_ok=True)\n",
    "\n",
    "# Read the input SDF file\n",
    "with open(input_file, 'r') as f:\n",
    "    data = f.read()\n",
    "\n",
    "# Split the SDF file into individual molecules\n",
    "molecules = data.split('$$$$\\n')\n",
    "\n",
    "# Iterate over each molecule and write it to a temporary SDF file\n",
    "for molecule in molecules:\n",
    "    # Extract the molecule title\n",
    "    title = molecule.split('\\n', 1)[0]\n",
    "\n",
    "    # Write the molecule to a temporary SDF file\n",
    "    with open('temp.sdf', 'w') as f:\n",
    "        f.write(molecule)\n",
    "\n",
    "    # Use obminimize to optimize the molecule geometry using MMFF94 force field\n",
    "    subprocess.run(['obminimize','-ff', 'MMFF94', '-c' '1e-6' '-n', '1000','temp.sdf'])\n",
    "\n",
    "    # Convert the optimized molecule to PDBQT format and write it to the output file\n",
    "    pdbqt_file = os.path.join(output_dir, f'{title}.pdbqt')\n",
    "    subprocess.run(['obabel', '-isdf', 'temp.sdf', '-opdbqt', '-O', pdbqt_file, '--partialcharge gasteiger'])\n",
    "\n",
    "# Remove the temporary SDF file\n",
    "os.remove('temp.sdf')\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "id": "974e763b-2de9-414c-8b42-4a4935e90ba6",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "1. ZINC000030663053 30663053 CSC025831865 MCULE-1695688201 MolPort-009-285-752 PubChem-135929469: -11.40\n",
      "2. Z358338076 ZINC000030663053 30663053 CSC025831865 MCULE-1695688201 MolPort-009-285-752 PubChem-135929469: -11.30\n"
     ]
    }
   ],
   "source": [
    "import os\n",
    "import subprocess\n",
    "\n",
    "# Define the input file and the Vina executable\n",
    "vina_exe = 'vina'\n",
    "\n",
    "# Define the input and output directories\n",
    "input_dir = 'Processed_Ligand'\n",
    "output_dir = 'Docked_Out'\n",
    "if not os.path.exists(output_dir):\n",
    "    os.makedirs(output_dir)\n",
    "\n",
    "# Define the number of top scoring molecules to keep\n",
    "top_n = 2\n",
    "\n",
    "# Run Vina for each ligand and save the log file with the molecule name\n",
    "log_files = []\n",
    "for ligand in os.listdir(input_dir):\n",
    "    if ligand.endswith('.pdbqt'):\n",
    "        # Run Vina for the ligand\n",
    "        ligand_path = os.path.join(input_dir, ligand)\n",
    "        ligand_name = ligand.split('.')[0]\n",
    "        output_file = os.path.join(output_dir, ligand_name + '_out.pdbqt')\n",
    "        log_file = os.path.join(output_dir, ligand_name + '_log.txt')\n",
    "        cmd = [vina_exe, '--config', 'config.txt', '--ligand', ligand_path,\n",
    "               '--out', output_file, '--log', log_file]\n",
    "        subprocess.run(cmd)\n",
    "        log_files.append(log_file)\n",
    "\n",
    "# Find the molecule with the lowest negative affinity\n",
    "affinities = []\n",
    "for log_file in log_files:\n",
    "    with open(log_file, 'r') as f:\n",
    "        # Find the line with the affinity value\n",
    "        for line in f:\n",
    "            if line.startswith('-----+------------+----------+----------'):\n",
    "                next(f)  # skip the separator row\n",
    "                affinity_line = next(f)  # the line with the affinity value\n",
    "                affinity = float(affinity_line.split()[1])\n",
    "                affinities.append((log_file, affinity))\n",
    "\n",
    "# Sort the affinities and print the top molecules\n",
    "affinities.sort(key=lambda x: x[1])\n",
    "top_molecules = []\n",
    "for i, (log_file, affinity) in enumerate(affinities[:top_n]):\n",
    "    molecule_name = os.path.basename(log_file).split('_log.txt')[0]\n",
    "    print(f'{i+1}. {molecule_name}: {affinity:.2f}')\n",
    "    top_molecules.append(f'{molecule_name}: {affinity:.2f}')\n",
    "    \n",
    "# Save the list of top molecules\n",
    "with open(os.path.join(output_dir, 'top_molecules.txt'), 'w') as f:\n",
    "    f.write('\\n'.join(top_molecules))\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "id": "0302443e-d8a7-49f4-a3ef-e0e88ab2fd13",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Converted Docked_Out\\Z358338076 ZINC000030663053 30663053 CSC025831865 MCULE-1695688201 MolPort-009-285-752 PubChem-135929469_out.pdbqt to Top_Docked_Molecules\\Z358338076 ZINC000030663053 30663053 CSC025831865 MCULE-1695688201 MolPort-009-285-752 PubChem-135929469.sdf\n",
      "Converted Docked_Out\\Z3243897341 ZINC000096124391 ZINC96124391 96124391 MCULE-8077027281 MolPort-030-077-218 PubChem-72719994_out.pdbqt to Top_Docked_Molecules\\Z3243897341 ZINC000096124391 ZINC96124391 96124391 MCULE-8077027281 MolPort-030-077-218 PubChem-72719994.sdf\n"
     ]
    }
   ],
   "source": [
    "import os\n",
    "import subprocess\n",
    "\n",
    "# Define the output directory\n",
    "output_dir = 'Top_Docked_Molecules'\n",
    "if not os.path.exists(output_dir):\n",
    "    os.makedirs(output_dir)\n",
    "\n",
    "# Read the top scoring compounds from the file\n",
    "top_molecules = []\n",
    "with open(r'Docked_Out\\top_molecules.txt', 'r') as f:\n",
    "    for line in f:\n",
    "        top_molecules.append(line.strip())\n",
    "\n",
    "# Convert the top scoring compounds to SDF format\n",
    "for molecule in top_molecules:\n",
    "    pdbqt_file = os.path.join('Docked_Out', molecule.split(':')[0] + '_out.pdbqt')\n",
    "    if not os.path.exists(pdbqt_file):\n",
    "        print(f'Warning: PDBQT file for {molecule} not found')\n",
    "        continue\n",
    "    sdf_file = os.path.join(output_dir, molecule.split(':')[0] + '.sdf')\n",
    "    cmd = f'obabel -i pdbqt {pdbqt_file} -o sdf -O {sdf_file}'\n",
    "    # Run the command and check the exit code to see if it succeeded\n",
    "    proc = subprocess.run(cmd, shell=True)\n",
    "    if proc.returncode == 0:\n",
    "        print(f'Converted {pdbqt_file} to {sdf_file}')\n",
    "    else:\n",
    "        print(f'Error converting {pdbqt_file} to SDF format')\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "f9116ac1-623c-4ec0-8772-55d1cdabbde4",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.12.7"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
