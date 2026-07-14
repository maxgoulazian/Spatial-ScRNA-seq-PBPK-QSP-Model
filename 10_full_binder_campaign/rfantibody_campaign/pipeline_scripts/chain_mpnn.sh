#!/bin/bash
# Wait for the LAST backbone campaign (costim_cea5) to appear, then finish, then run ProteinMPNN on all 4 arms.
LOG=/media/balthasar-lab/RAID1/RFantibody/campaign_out/chain.log
# poll up to 24h for the cea5 container to be created
for i in $(seq 1 1440); do docker ps -a --format '{{.Names}}' | grep -q '^costim_cea5$' && break; sleep 60; done
docker wait costim_cea5 >/dev/null 2>&1
sleep 5
echo "[$(date +%H:%M:%S)] all backbones done -> launching ProteinMPNN (all 4 arms)" >> $LOG
docker run -d --rm --gpus all --name costim_mpnn \
  -v /media/balthasar-lab/RAID1/RFantibody:/home \
  --entrypoint /bin/bash ullahsamee/rfantibody:latest \
  -lc "bash /home/run_proteinmpnn.sh > /home/campaign_out/mpnn.log 2>&1" >/dev/null 2>&1
echo "[$(date +%H:%M:%S)] costim_mpnn container launched" >> $LOG
