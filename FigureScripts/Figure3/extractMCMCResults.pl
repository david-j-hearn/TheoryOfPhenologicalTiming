#!/usr/bin/perl
#

@files = (
	"screenOutput.txt"
	);


$sep = 10;
foreach $file (@files)
	{
	open(IN, "$file");
	open(OUT, ">$file.Cleaned.ForTracer.txt");
	($species, $j, $part, $j) = split(/\./,$file);

	print OUT "gen\tonset.shape1\tonset.shape2\tduration.shape1\tduration.shape2\tlk\tprior\tpost\tmean.onset\tvariance.onset\tmean.duration\tvariance.duration\n";
	#print OUT "Gen\tLnPr\tLnL\tLnPosterior\tOnset.Shape1\tOnset.Shape2\tDuration.Shape1\tDuration.Shape2\tOnset.Acceptance\tDuration.Acceptance\tOnset.Mean\tOnset.Variance\tDuration.Mean\tDuration.Variance\n";

	$itP = 0;
	while($line = <IN>)
		{
		chomp($line);

		#[1] "i, (onset), (duration), accept, pr, lk, post, oRate, dRate  1 ( 56.882924683244 134.460151880806 ) ( 82.5096231888601 1325.71502221586 ) 0 0 0 0 -35.2149448778513 -1717.48141439572 -1752.69635927357 0.0025 0.0025 o mean 108.508067719054 o var 144.696884621496 d mean 15.0281794545652 d var 2.57500507181769"
		
		($j, $j, $j, $j, $j, $j, $j, $j, $j, $j, $it, $j, $oShape1, $oShape2, $j, $j, $dShape1, $dShape2, $j, $aO, $j, $aD, $j, $pr, $lk, $post, $j, $j, $j, $j, $oMean, $j, $j, $oVar, $j, $j, $dMean, $j, $j, $dVar) = split(/\s+/, $line); 

		if($dVar) 
			{
				#print OUT "$species\t$part\t$it\t$pr\t$lk\t$post\t$oShape1\t$oShape2\t$dShape1\t$dShape2\t$aO\t$aD\n";
			$dVar =~ s/"//eg;
			if($it - $itP > $sep)
				{
				$tit = $itP+$sep;
				while($tit < $it)
					{
					print OUT "$tit\t$oShape1\t$oShape2\t$dShape1\t$dShape2\t$lk\t$pr\t$post\t$oMean\t$oVar\t$dMean\t$dVar\n";
					#print OUT "$tit\t$pr\t$lk\t$post\t$oShape1\t$oShape2\t$dShape1\t$dShape2\t$aO\t$aD\t$oMean\t$oVar\t$dMean\t$dVar\n";
					$tit = $tit + $sep;
					}
				}
			print OUT "$it\t$oShape1\t$oShape2\t$dShape1\t$dShape2\t$lk\t$pr\t$post\t$oMean\t$oVar\t$dMean\t$dVar\n";
			#print OUT "$it\t$pr\t$lk\t$post\t$oShape1\t$oShape2\t$dShape1\t$dShape2\t$aO\t$aD\t$oMean\t$oVar\t$dMean\t$dVar\n";
			$itP = $it;
			}
		}
	close(IN);
	close(OUT);
	}


