create new chothia numbered input:
python pdb_to_chothia.py --pdb /media/balthasar-lab/RAID1/RFantibody/input/EpAb.pdb --csv /media/balthasar-lab/RAID1/RFantibody/my_new_seqs.csv --out /media/balthasar-lab/RAID1/RFantibody/input/EpAb_chothia.pdb


Start RFantibody container:
docker run --gpus all -v .:/home --memory 600g -it ullahsamee/rfantibody

Input Prep:
(inside the container)
poetry run python /home/scripts/util/chothia2HLT.py -H H -L L -T C '/home/input/8bw0_chothia.pdb'

Run: (edit antibody_pdbdesign.sh with specifics)
bash /home/scripts/examples/rfdiffusion/nanobody_pdbdesign.sh 
