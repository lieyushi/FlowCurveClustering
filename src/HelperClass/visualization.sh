Rscript correlation_sheet.r gamma_ranking 1_colnames gamma_ranking.pdf

Rscript correlation_sheet.r dbindex_ranking 2_colnames dbindex_ranking.pdf

Rscript correlation_sheet.r validity_ranking 3_colnames validity_ranking.pdf

Rscript correlation_sheet.r silhouette_ranking 0_colnames silhouette_ranking.pdf

Rscript highlighting_best.r gamma_best 1_colnames gamma_best.pdf

Rscript highlighting_best.r dbindex_best 2_colnames dbindex_best.pdf

Rscript highlighting_best.r validity_best 3_colnames validity_best.pdf

Rscript highlighting_best.r silhouette_best 0_colnames silhouette_best.pdf

convert gamma_ranking.pdf gamma_ranking.png

convert silhouette_ranking.pdf silhouette_ranking.png

convert validity_ranking.pdf validity_ranking.png

convert dbindex_ranking.pdf dbindex_ranking.png

convert -transparent white gamma_best.pdf gamma_best.png

convert -transparent white dbindex_best.pdf dbindex_best.png

convert -transparent white validity_best.pdf validity_best.png

convert -transparent white silhouette_best.pdf silhouette_best.png

if [ -f gamma_std ]; then

	Rscript highlighting_best.r gamma_std 1_colnames gamma_std.pdf

	Rscript highlighting_best.r dbindex_std 2_colnames dbindex_std.pdf

	Rscript highlighting_best.r validity_std 3_colnames validity_std.pdf

	Rscript highlighting_best.r silhouette_std 0_colnames silhouette_std.pdf
	
	mv *_std $1
	mv *_std.pdf $1
fi

mv *_colnames $1
mv *_ranking.pdf $1
mv *_ranking $1
mv assembled $1
mv *_best.pdf $1
mv *.png $1
mv *_best $1
