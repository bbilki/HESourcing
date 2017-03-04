#!/bin/sh

if [ ! -d Plots ]
then
	mkdir Plots
fi

if [ ! -d Plots/$1 ]
then
	mkdir Plots/$1
fi

if [ ! -d Histos ]
then
	mkdir Histos
fi

if [ ! -d Histos/$1 ]
then
	mkdir Histos/$1
fi

if [ ! -d NTuples ]
then
	mkdir NTuples
fi

if [ ! -d NTuples2 ]
then
	mkdir NTuples2
fi

cmsRun hesourcing_cfg.py $1 2
wait
mv N_$1.root NTuples
wait
rm Data/USC_$1.root
 
cd Files/reNTupler
./reNTupler $1
wait
mv N2_$1.root ../../NTuples2
mv H2_$1.root ../../Histos/$1
wait
cd -

cd Files
./plotter $1
wait
mv *.png ../Plots/$1/
wait
cd -
# 
# # cd Plots/$1
# # convert -delay 100 -loop 0 PH*.png PH.gif
# # wait
# # cd -
# 
cd Plots
tar -cf plots.tar $1
wait
gzip plots.tar
wait
scp -P 53222 plots.tar.gz bbilki@feynman.physics.uiowa.edu:/var/www/html/HESourcing
wait
# ssh -p 53222 bbilki@feynman.physics.uiowa.edu "cd /var/www/html/HESourcing ; tar -zxf /var/www/html/HESourcing/plots.tar.gz ; rm plots.tar.gz ; ./makePlotList.sh"
ssh -p 53222 hfSX5@feynman.physics.uiowa.edu "cd /var/www/html/HESourcing ; tar -zxf /var/www/html/HESourcing/plots.tar.gz ; rm plots.tar.gz ; ./makePlotList.sh"
wait
rm plots.tar.gz
cd -

# rm Files/MeanADCperTile.txt
# cd Files
# for i in 286975 286977 286979 286981 286983 286985 286987 286989 287083 287084 287085 287086 287087
# do
# 	./plotter $i
# 	wait
# 	mkdir ../Plots/$i
# 	wait
# 	mv *.png ../Plots/$i
# 	wait
# 	cd ../Plots/$i
# 	convert -delay 100 -loop 0 PH*.png PH.gif
# 	wait
# 	cd ../../Files
# done
# 
# cd ../Plots
# # cd Plots
# tar -cf plots.tar *
# wait
# gzip plots.tar
# wait
# scp -P 53222 plots.tar.gz bbilki@feynman.physics.uiowa.edu:/var/www/html/HESourcing
# wait
# ssh -p 53222 bbilki@feynman.physics.uiowa.edu "cd /var/www/html/HESourcing ; tar -zxf /var/www/html/HESourcing/plots.tar.gz ; rm plots.tar.gz ; ./makePlotList.sh"
# wait
# rm plots.tar.gz
# cd -




