#!/bin/bash
# wait for the 41BB->CD27 campaign container to finish
docker wait costim_campaign >/dev/null 2>&1
sleep 5
mkdir -p /media/balthasar-lab/RAID1/RFantibody/campaign_out/CD3
docker run -d --rm --gpus all --name costim_cd3 \
  -v /media/balthasar-lab/RAID1/RFantibody:/home \
  --entrypoint /bin/bash ullahsamee/rfantibody:latest \
  -lc "bash /home/run_cd3_campaign.sh" >/dev/null 2>&1
echo "[$(date +%H:%M:%S)] CD3 container launched after CD27 completion" >> /media/balthasar-lab/RAID1/RFantibody/campaign_out/chain.log
