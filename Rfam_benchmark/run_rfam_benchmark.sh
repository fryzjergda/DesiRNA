#!/bin/bash

#for d in ./*/; do 
#	cd "$d" && for i in *txt; do DesiRNA.py -f $i -t 2 ;done
#	cd ../	
#done

tim=60

(cd rfam_01-05 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > rfam_01-05_o.bench && cat *$tim'_new.best_str' > rfam_01-05_new.bench &&    cp *.bench ../.) &
(cd rfam_06-10 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > rfam_06-10_o.bench && cat *$tim'_new.best_str' > rfam_06-10_new.bench &&    cp *.bench ../.) &
(cd rfam_11-15 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > rfam_11-15_o.bench && cat *$tim'_new.best_str' > rfam_11-15_new.bench &&    cp *.bench ../.) &
(cd rfam_16-20 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > rfam_16-20_o.bench && cat *$tim'_new.best_str' > rfam_16-20_new.bench &&    cp *.bench ../.) &
(cd rfam_21-25 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > rfam_21-25_o.bench && cat *$tim'_new.best_str' > rfam_21-25_new.bench &&    cp *.bench ../.) &
(cd rfam_26-30 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > rfam_26-30_o.bench && cat *$tim'_new.best_str' > rfam_26-30_new.bench &&    cp *.bench ../.) &



