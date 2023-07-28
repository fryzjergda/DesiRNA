#!/bin/bash

#for d in ./*/; do 
#	cd "$d" && for i in *txt; do DesiRNA.py -f $i -t 2 ;done
#	cd ../	
#done

tim=600
sf=dmt


(cd ete_V1_01-10 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim -acgu on -p 1999 -R 10 -l 1000 -e 10 -sf $sf >&/dev/null;done && cat *$tim.best_str > ete_V1_1-10_o.bench && cat *$tim'_new.best_str' > ete_V1_1-10_new.bench &&    cp *.bench ../.) &
(cd ete_V1_11-20 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim -acgu on -p 1999 -R 10 -l 1000 -e 10 -sf $sf >&/dev/null;done && cat *$tim.best_str > ete_V1_11-20_o.bench && cat *$tim'_new.best_str' > ete_V1_11-20_new.bench &&    cp *.bench ../.) &



