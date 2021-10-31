#!/bin/bash
par=".par"
pdf=".pdf"
txt=".txt"
post=".post"
res=".results"
co=".corner"
ch=".chain"
chp="chp_"
dbr="dbr_"
trh="trh_"
all="_eall"
for timfile in chp_*.tim
do
    name="$(cut -d'.' -f1 <<<"$timfile")"
    pulsar="$(cut -d'_' -f2 <<<"$name")"
    echo "$pulsar"
    /local/scratch/yangliu/glitch/run_enterprise/run_enterprise.py --gl-all --auto-add --emcee -N 4000 --nwalkers 128 -t 2 "$dbr$pulsar$par" "$timfile" --truth-file "$trh$pulsar$txt" --plot-chain --plot-derived -j --red-prior-log  -A -8 --tspan-mult 1.1 --glitch-alt-f0 --glitch-alt-f0t 200 --alt-f0t-gltd --glitch-epoch-range 100 --measured-prior --measured-sigma 50 --glitch-td-min 1 --glitch-td-max 3 --glitch-f0d-range 3.0 --glitch-f0-range 0.8 --glitch-f1-range 0.8 --glitch-td-split 2
    mv "$pulsar$co$pdf" "$pulsar$all$co$pdf"
    mv "$pulsar$ch$pdf" "$pulsar$all$ch$pdf"
    mv "$dbr$pulsar$par$post" "$dbr$pulsar$all$par$post"
    mv "$dbr$pulsar$par$res" "$dbr$pulsar$all$par$res"
done
