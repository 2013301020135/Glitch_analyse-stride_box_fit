#!/bin/bash
par=".par"
pdf=".pdf"
txt=".txt"
post=".post"
res=".results"
co=".corner"
tr=".traceplot"
ru=".runplot"
chp="chp_"
dbr="dbr_"
trh="trh_"
all="_dall"
for timfile in chp_*.tim
do
    name="$(cut -d'.' -f1 <<<"$timfile")"
    pulsar="$(cut -d'_' -f2 <<<"$name")"
    echo "$pulsar"
    /local/scratch/yangliu/glitch/run_enterprise/run_enterprise.py --gl-all --auto-add --dynesty --nlive 500 -t 2 "$dbr$pulsar$par" "$timfile" --truth-file "$trh$pulsar$txt" --dynesty-plots --plot-derived -j --red-prior-log  -A -8 --tspan-mult 1.1 --glitch-alt-f0 --glitch-alt-f0t 200 --alt-f0t-gltd --glitch-epoch-range 100 --measured-prior --measured-sigma 50 --glitch-td-min 1 --glitch-td-max 3 --glitch-f0d-range 3.0 --glitch-f0-range 0.8 --glitch-f1-range 0.8 --glitch-td-split 2
    mv "$pulsar$co$pdf" "$pulsar$all$co$pdf"
    mv "$pulsar$tr$pdf" "$pulsar$all$tr$pdf"
    mv "$pulsar$ru$pdf" "$pulsar$all$ru$pdf"
    mv "$dbr$pulsar$par$post" "$dbr$pulsar$all$par$post"
    mv "$dbr$pulsar$par$res" "$dbr$pulsar$all$par$res"
done
