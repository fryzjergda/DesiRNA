#!/bin/bash

#for d in ./*/; do 
#	cd "$d" && for i in *txt; do DesiRNA.py -f $i -t 2 ;done
#	cd ../	
#done

tim=60

(cd ete_V2_01-10 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > ete_V2_1-10_o.bench && cat *$tim'_new.best_str' > ete_V2_1-10_new.bench &&    cp *.bench ../.) &
(cd ete_V2_11-20 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > ete_V2_11-20_o.bench && cat *$tim'_new.best_str' > ete_V2_11-20_new.bench &&    cp *.bench ../.) &
(cd ete_V2_21-30 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > ete_V2_21-30_o.bench && cat *$tim'_new.best_str' > ete_V2_21-30_new.bench &&    cp *.bench ../.) &
(cd ete_V2_31-40 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > ete_V2_31-40_o.bench && cat *$tim'_new.best_str' > ete_V2_31-40_new.bench &&    cp *.bench ../.) &
(cd ete_V2_41-50 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > ete_V2_41-50_o.bench && cat *$tim'_new.best_str' > ete_V2_41-50_new.bench &&    cp *.bench ../.) &
(cd ete_V2_51-60 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > ete_V2_51-60_o.bench && cat *$tim'_new.best_str' > ete_V2_51-60_new.bench &&    cp *.bench ../.) &
(cd ete_V2_61-70 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > ete_V2_61-70_o.bench && cat *$tim'_new.best_str' > ete_V2_61-70_new.bench &&    cp *.bench ../.) &
(cd ete_V2_71-80 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > ete_V2_71-80_o.bench && cat *$tim'_new.best_str' > ete_V2_71-80_new.bench &&    cp *.bench ../.) &
(cd ete_V2_81-90 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > ete_V2_81-90_o.bench && cat *$tim'_new.best_str' > ete_V2_81-90_new.bench &&    cp *.bench ../.) &
(cd ete_V2_91-100 &&  for i in *txt; do nohup DesiRNA.py -f $i -t $tim >&/dev/null;done && cat *$tim.best_str > ete_V2_91-100_o.bench && cat *$tim'_new.best_str' > ete_V2_91-100_new.bench &&    cp *.bench ../.) &



