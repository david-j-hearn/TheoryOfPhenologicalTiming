#!/usr/bin/perl

#must install the following Perl module to use this script
use Proc::Background;

@cmds = (); 
$nProc = 26; #set to the number of processors available on your Linux machine


@files = ("AnemoneQuinquefolia.txt", "DicentraCucullaria.txt", "PrimulaMeadia.txt", "CamassiaScilloides.txt", "EnemionBiternatum.txt", "SanguinariaCanadensis.txt", "CardamineConcatenata.txt","ErythroniumAmericanum.txt","ThalictrumThalictroides.txt", "ClaytoniaVirginica.txt", "MertensiaVirginica.txt", "CollinsiaVerna.txt", "PodophyllumPeltatum.txt");



if(scalar @cmds ==0 )
{
	mkdir("../MCMCReplicates");
	foreach $file (@files)
	{
		$cmd = "Rscript --no-restore --no-save -e 'source(\"./runModelFit_SplitCentury_MCMC.R\"); run_BayesianReplicate(file=\"$file\", firstHalf=T)' &> ../MCMCReplicates/$file.FirstHalf.SO.redo.txt";
		push(@cmds, $cmd);
		$cmd = "Rscript --no-restore --no-save -e 'source(\"./runModelFit_SplitCentury_MCMC.R\"); run_BayesianReplicate(file=\"$file\", firstHalf=F)' &> ../MCMCReplicates/$file.SecondHalf.SO.redo.txt";
		push(@cmds, $cmd);
	}
}


$totCmds = scalar @cmds;
print "There are " . $totCmds . " commands to run in total.\n";

for($i = 0; $i < $nProc; $i++)
{
	$cmd = pop(@cmds);
	print "Running command " . $cmd . "\n";
	$proc[$i] = Proc::Background->new("$cmd");
}


$totCnt = $nProc;
while( scalar @cmds > 0 )
{
	for($i = 0; $i< $nProc; $i++)
	{
		if(!$proc[$i]->alive)
		{
			$cmd = pop(@cmds);
			print "Running command " . $cmd . " of $totCmds with " . scalar @cmds . " commands left.\n";
			$proc[$i] = Proc::Background->new("$cmd");
		}
		if(scalar @cmds == 0 ) 
		{
			print "Out of commands\n";
			$done = 0;
			while( $done < $nProc )
			{
				$done = 0;
				for($j = 0; $j< $nProc; $j++)
				{
					if(!$proc[$j]->alive)
					{
						$done++;
					}
				}
				if($done < $nProc) 
				{
					sleep(10);
				}
			}
			exit(0);
		}
	}
	sleep(10);
}

if(scalar @cmds == 0 ) 
{
	print "Out of commands\n";
	$done = 0;
	while( $done < $nProc )
	{
		$done = 0;
		for($j = 0; $j< $nProc; $j++)
		{
			if(!$proc[$j]->alive)
			{
				$done++;
			}
		}
		if($done < $nProc) 
		{
			sleep(10);
		}
	}
	exit(0);
}
