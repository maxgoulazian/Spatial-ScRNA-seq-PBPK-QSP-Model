#!/bin/bash
# wait for CD3 container; if not yet created, poll for it first
for i in $(seq 1 720); do docker ps -a --format '{{.Names}}' | grep -q '^costim_cd3$' && break; sleep 60; done
docker wait costim_cd3 >/dev/null 2>&1
sleep 5
mkdir -p /media/balthasar-lab/RAID1/RFantibody/campaign_out/CEA5
docker run -d --rm --gpus all --name costim_cea5 \
  -v /media/balthasar-lab/RAID1/RFantibody:/home \
  --entrypoint /bin/bash ullahsamee/rfantibody:latest \
  -lc "bash /home/run_cea5_campaign.sh" >/dev/null 2>&1
echo "[$(date +%H:%M:%S)] CEA5 container launched after CD3 completion" >> /media/balthasar-lab/RAID1/RFantibody/campaign_out/chain.log
