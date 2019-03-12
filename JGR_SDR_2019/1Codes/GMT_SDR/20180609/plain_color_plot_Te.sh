filename=global_Te_SDR_20180609_Te_linear.ps
gmt pscoast -R-169/190/-58/78 -JM6i -P -Ba -Ggrey -K > $filename
gmt psxy data.txt -R -J -P -Sc.2 -CGMT_seisTe.cpt -K -O >> $filename
gmt psxy data_eg.txt -R -J -P -Sa.6 -CGMT_seisTe.cpt -K -O >> $filename
gmt psxy data_EAR.txt -R -J -P -Si.2 -CGMT_seisTe.cpt -K -O >> $filename
gmt psscale -Dx0i/5.8i+w6i/0.2i+h -CGMT_seisTe.cpt -B+tcontinuous -O >> $filename
gmt psconvert -Tf $filename
