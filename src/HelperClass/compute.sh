./main $1
mv *_colnames R\ implementation/
mv *_ranking R\ implementation/
mv *_best R\ implementation/
mv assembled R\ implementation/
if [ -f gamma_std ]
then
	mv *_std R\ implementation/
fi
