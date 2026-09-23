#!/bin/bash
HOSTIP=$(ip route show default | awk '{print $3}')
echo "host=$HOSTIP"
for port in 10808 10809 7890; do
  timeout 15 curl -s -x "http://$HOSTIP:$port" -r 0-3000000 -o /dev/null \
    -w "proxy $port: %{speed_download} B/s http=%{http_code}\n" \
    https://hgdownload.soe.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw
done
echo "--- direct (no proxy) 5s sample ---"
timeout 6 curl -s -r 0-3000000 -o /dev/null \
  -w "direct: %{speed_download} B/s http=%{http_code}\n" \
  https://hgdownload.soe.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw
