#!/usr/bin/perl




#12Fev2013
# - allow any unit (species, station, ...) superior or equal to the used unit level
#   to be selected using the -selectedunits options:
# 		-unit indiv -selectedunits station_file
# 	with station_file looking like that:
# 	sp1|LOCATIONX
# 	sp2|LOCATIONY
# - use File::Basename to get file name and path
#18Avril 2012

#10Jul2013
# - add the RNAseq paradigm

#26Sept2013
# - reverse the noW option: now becomes default and is replaced by -useW

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use YAML;
use Pod::Usage;
use File::Basename;

####TODO
# - fix the W that produces errors (2 seq identical but one W)
# - fix the NCBI format names
# - allow some flexibility for the unit names so that it fits into a strict name... (or change the concatenate script)



##### input
# gene files, the following convention: NAME_date.[fas|mase]
# there is two categories of gene files: master genes or optional genes
# what about three categories of genes: must be, as much as possible (with a minimum ratio), and optional? (may be with 1, 2 and 3 flags) --> not for the moment


#### options
# by definition the optional genes can be present or not
# the master genes must be present at a given ratio (eg 2 out of 3)

# the assembly can be done at 2 unit levels: individuals or station (within a species)

# by default, only one seq per unit is kept

# if several sequence are available, the seq that has the maximum [AGCT] is selected
# sequences with more than a given ambiguity content are rejected

####output
# one fasta sequence per gene with the selected sequence
# an assembly table, that will be eventually used to produce concatenation

my @masterGenes;
my @optionalGenes;
my $minOptionalGenes = 0;
my $unit = 'indiv';
my $outtable = 'assembly_table.tsv';
my $prefix = 'chimera_assembler_';
my $maxAmbiguity = 5;
my $noW = 1;
my $useW = 0;
my $outTab = 'assembly_table.tsv';
my $criteria = 'bestindiv';
my $selectedUnitsFile;
my $exportseq;
my $motuFile;
my $help;
my $NoTranscript;
my $NoNCBI;



GetOptions(	"master=s" 		=> \@masterGenes,
		"optional=s" 		=> \@optionalGenes,
		"min=i" 			=> \$minOptionalGenes,
		"unit=s" 			=> \$unit,
		"prefix=s" 			=> \$prefix,
		"maxambiguity=f" 	=> \$maxAmbiguity,
		"useW" 				=> \$useW,
		"outtab=s" 			=> \$outTab,
		"criteria=s" 		=> \$criteria,
		"selectedunits=s" 	=> \$selectedUnitsFile,
		"exportseq" 		=> \$exportseq,
		"motu=s" 			=> \$motuFile,
		"help" 				=> \$help,
		"noTranscript" 		=> \$NoTranscript,
		"noNCBI" 			=> \$NoNCBI,
);

pod2usage(-exitstatus => 0, -verbose => 2) if $help;

if(@masterGenes == 0 && @optionalGenes == 0) {
# 	print $usage;
	pod2usage(-exitstatus => 0, -verbose => 1);
}

if(!$unit =~ /^(station|indiv|species|motu)$/) {
	pod2usage(-exitstatus => 1, -verbose => 1);
}

if(!$criteria =~ /^(bestindiv|bestchimera)$/) {
	pod2usage(-exitstatus => 1, -verbose => 1);
}

if($unit eq 'motu' and !$motuFile) {
	pod2usage(-exitstatus => 1, -verbose => 1);
}

$noW = 0 if($useW);


open OUT, ">ind_hash.yaml";
open OUT2, ">yas_hash.yaml";
open OUT3,  ">motu2seq.yaml";

#make a total list of genes
my @allGenesFiles = (@masterGenes, @optionalGenes);
my @allGenes;
# foreach (@allGenesFiles) {
# 	s/_\d+\.\w+//;
# # 	print "gene : $_\n";
# 	push @allGenes, $_;
# }

my %geneCategories;
my %geneToGeneFile;
foreach my $file (@masterGenes) {
	 my ($name,$path) = fileparse($file);
	(my $gene = $name) =~ s/_\d+\.\w+$//;
	push @allGenes, $gene;
	$geneCategories{$gene} = 'master';
	$geneToGeneFile{$gene} = $file;
}
foreach my $file (@optionalGenes) {
	my ($name,$path) = fileparse($file);
	(my $gene = $name) =~ s/_\d+\.\w+$//;
	push @allGenes, $gene;
	$geneCategories{$gene} = 'optional';
	$geneToGeneFile{$gene} = $file;
}
my $nRequiredMaster = scalar @masterGenes;

# print join("\n", @allGenesFiles), "\n";

#####List of units to keep
my %unitToKeep;
if($selectedUnitsFile) {
	open LIST, $selectedUnitsFile;
	while (my $i = <LIST>) {
		chomp $i;
		#not a fan of this hack, but well... here I escape the | so that I can do some grep on it
		$i =~ s:\|:\\|:;
		$i =~ s:\[:\\[:;
		$i =~ s:\]:\\]:;
		$unitToKeep{$i} = 1;
	}
}


############# HASH the seq
my $indHash;
my $yas2info;


foreach my $genefile (@allGenesFiles) {
# 	(my $gene = $genefile) =~ s/_\d+\.\w+//;
	my ($name,$path) = fileparse($genefile);
	(my $gene = $name) =~ s/_\d+\.\w+$//;
	my $seqin = Bio::SeqIO->new( -file => $genefile, -format => 'fasta', -alphabet => 'dna');
	my $n = 0;
	while(my $seq = $seqin->next_seq ) {
		my $id = $seq->id;
# 		print "$gene $id\n";
		last if $id =~ /^[_-]+Primer/i;
		last if $id =~ /^[_-]+RESERVE/i;
		next if $id =~ /^--/;
		next if $id  =~ /^__/;
		next if ($noW && $id =~ /^W_/);
		unless($id =~ /\|\w+\|/) {
			print "Don't get this one: $id\n";
			next;
		}
		++$n;

		my $seqinfo = parse_formatv4($id);
		my $ind = $seqinfo->{ind};
		my $station = $seqinfo->{loc};
		my $date = $seqinfo->{date};
		my $sp = $seqinfo->{head};
		my $firstyas = $seqinfo->{firstyas};
		my $seqType = $seqinfo->{seqType};

		next if($NoTranscript and $seqType eq 'transcript');
		next if($NoNCBI and $seqType eq 'NCBI');


# 		print "$sp $station $date $ind\n";
		#remove the working flag
		$sp =~ s/^W_//;
		my @YASes = @{$seqinfo->{yas}};
		if(exists $indHash->{$sp}->{$station}->{$date}->{$ind}->{$gene}) {
			print "Pas content avec le $gene: la meme cellule a plusieur prétendants:
			$indHash->{$sp}->{$station}->{$date}->{$ind}->{$gene}->{id}
			$id\n";
			next;
		}

		#initialize the gene count if not yet done
		if (!defined $indHash->{$sp}->{$station}->{$date}->{$ind}->{nmaster}) {
			$indHash->{$sp}->{$station}->{$date}->{$ind}->{nmaster} = 0;
		}
		if (!defined $indHash->{$sp}->{$station}->{$date}->{$ind}->{noptional}) {
			$indHash->{$sp}->{$station}->{$date}->{$ind}->{noptional} = 0;
		}

		#get the sequence quality parameters
		my $seqstring = $seq->seq();
		my $seqscore = calculate_seqscore($seqstring);
		my $ambiguityscore = calculate_ambiguityscore($seqstring);
		my $ambiguityRate = 100 * $ambiguityscore / $seqscore;

		#skip this if not good quality
		next if($ambiguityRate > $maxAmbiguity);

		#fill the database
		$indHash->{$sp}->{$station}->{$date}->{$ind}->{$gene}->{id} = $id;
		push @{$indHash->{$sp}->{$station}->{$date}->{$ind}->{$gene}->{yas}}, @YASes;
		$indHash->{$sp}->{$station}->{$date}->{$ind}->{$gene}->{'seqscore'} = $seqscore;
		$indHash->{$sp}->{$station}->{$date}->{$ind}->{$gene}->{ambiguityscore} = $ambiguityRate;

		#fill the reverse hash (to go from the yas to the seq info)
		if(defined $yas2info->{$firstyas}->{$ind}->{$date}) {
			print "WTF: the same YAS encountered several times for the same individual! $firstyas\n";
		}
		else { #have to go pretty deep due to the same YAS used several times for different ind...
			$yas2info->{$firstyas}->{$ind}->{$date} = {
				sp => $sp,
				loc => $station,
				gene => $gene,
			};
		}

		#increment the gene count for this ind
		if($geneCategories{$gene} eq 'master') {
			++$indHash->{$sp}->{$station}->{$date}->{$ind}->{nmaster};
		}
		elsif($geneCategories{$gene} eq 'optional') {
			++$indHash->{$sp}->{$station}->{$date}->{$ind}->{noptional};
		}
		else {
			print "I don't know to which category this gene $gene should be assigned\n";
		}
	}
}

####### Get the motu file if available
my %motu2seq;
if($motuFile) {
	open MOTU, $motuFile;
	while (<MOTU>) {
		chomp;
		my @f = split "\t";
		my $seq = $f[0];
		my $m = $f[1];
		push @{$motu2seq{$m}}, $seq;
		my ($sp, $ind, $loc, $date, $yas);
		#test wether this is a complete name or a core name
		if($seq =~ /\|.*\|/) {
			my $desc = parse_formatv4($seq);
			$ind = $desc->{ind};
			$loc = $desc->{loc};
			$date = $desc->{date};
			$sp = $desc->{sp};
			$yas = $desc->{firstyas};
		}
		elsif($seq =~ /(\w+)_(\d+)_(\w+)_(\w+)/) {
			$loc = $1;
			$date = $2;
			$ind = $3;
			$yas = $4;
			#carefull with the NCBI seq (apply the AN to the ind)
#			if($ind eq 'NCBI' or $ind eq 'RNAseq') {
#				$ind = $yas;
#			}
		}
		else {
			print "I don't get the following MOTU seq format: $seq\n";
			next;
		}

		#from the yas, find the corresponding individual, and assign the MOTU
		if(exists $yas2info->{$yas}->{$ind}->{$date}) {
			 my $sp_y = $yas2info->{$yas}->{$ind}->{$date}->{sp};
			 #fast checking that we are speaking about the same ind_hash
			my $loc_y = $yas2info->{$yas}->{$ind}->{$date}->{loc};
			if($loc_y ne $loc) {
				print "WTF: the MOTU seq YAS ($seq) points to an individual that is incongruent with the original seq ($loc_y vs $loc)\n";
				next;
			}
			#add the MOTU info to the big hash
			if(defined $indHash->{$sp_y}->{$loc}->{$date}->{$ind}->{MOTU}) {
				if($indHash->{$sp_y}->{$loc}->{$date}->{$ind}->{MOTU} eq $m) {next }
				else {
					print "WTF: already defined and different MOTU for this ind: $seq\n";
					next;
				}
			}
			else {
				$indHash->{$sp_y}->{$loc}->{$date}->{$ind}->{MOTU} = $m;
			}
		}
		else {
			print "Can't connect this MOTU seq to any known individual: $seq\n";
			next;
		}

	}
}

#### afterward find the pop with a single MOTU -> apply it to the other individuals
if($motuFile) {

	foreach my $sp (sort keys %$indHash) {
		foreach my $station (sort keys %{$indHash->{$sp}}) {
			my %nmotu;
			foreach my $date (sort keys %{$indHash->{$sp}->{$station}}) {
				foreach my $ind (sort keys %{$indHash->{$sp}->{$station}->{$date}}) {
					if(defined $indHash->{$sp}->{$station}->{$date}->{$ind}->{MOTU}) {
						my $m = $indHash->{$sp}->{$station}->{$date}->{$ind}->{MOTU};
						$nmotu{$m} = 1;
					}
				}
			}
			#how many motus in this pop?
			my @M = keys %nmotu;
			my $n = scalar @M;
			if($n == 0) { next }
			#if a single MOTU, apply it to the individuals
			elsif($n == 1) {
				my $Motu = $M[0];
				foreach my $date (sort keys %{$indHash->{$sp}->{$station}}) {
					foreach my $ind (sort keys %{$indHash->{$sp}->{$station}->{$date}}) {
						if(!defined $indHash->{$sp}->{$station}->{$date}->{$ind}->{MOTU}) {
							$indHash->{$sp}->{$station}->{$date}->{$ind}->{MOTU} = $Motu;
# 							print "$sp $station $date $ind <- $Motu\n";
						}
					}
				}
			}
			else {
				print "Can't propagate MOTU inside this pop: $sp $station (", join(',', @M), ")\n";
			}
		}
	}
}


print OUT Dump($indHash);
print OUT2 Dump($yas2info);
print OUT3 Dump(%motu2seq);





open TAB, ">$outTab";

print TAB "unit_name\t", join("\t", @allGenes) , "\n";

##############making the assembly table


my $nIncluded = 0;
my $nRejected = 0;
my %seqToExport;

###unit is individual
#just print out all the seq that fit the master and optinal gene cutoffs

if($unit eq 'indiv') {
	#find the individuals that fit the minimum requirements
	foreach my $sp (sort keys %$indHash) {
		foreach my $station (sort keys %{$indHash->{$sp}}) {
			foreach my $date (sort keys %{$indHash->{$sp}->{$station}}) {
				foreach my $ind (sort keys %{$indHash->{$sp}->{$station}->{$date}}) {
					my $unitId = "$sp|${station}_${date}_${ind}";
					#skip if a non selected unit
					if($selectedUnitsFile) {
						unless(map($unitId =~ /$_/, keys %unitToKeep)) { next }
					}
					my $localHash =  \%{$indHash->{$sp}->{$station}->{$date}->{$ind}};
					if($localHash->{nmaster} == $nRequiredMaster && $localHash->{noptional} >= $minOptionalGenes ) {
						print TAB "$sp|${station}_${date}_$ind";
						++$nIncluded;
# 						print "including this ind: $sp $station $date $ind\n";
						foreach my $gene (@allGenes) {
							#####
# 							print "$gene\n";
							if(exists $localHash->{$gene}) {
								my $seqid = $localHash->{$gene}->{id};
								print TAB "\t$seqid";
								$seqToExport{$gene}{$seqid} = 1;
							}
							else {
								print TAB "\t";
							}
						}
						print TAB "\n";
					}
					else {
						++$nRejected;
					}
				}
			}
		}
	}
}

###The unit is the species-station
#either select the best indivial per station, or make the best chimera

elsif($unit eq 'station') {
	foreach my $sp (sort keys %$indHash) {
		foreach my $station (sort keys %{$indHash->{$sp}}) {
			#this is the unit level
			my $unitId = "$sp|${station}";
			#skip if a non selected unit
			if($selectedUnitsFile) {
				unless(map($unitId =~ /$_/, keys %unitToKeep)) { next }
			}
			#get all the individual below this level, and either select the best indiv or the best chimera
			if($criteria eq 'bestindiv') {
				my $bestIndHash;
				my $bestACGTScore = 0;
				foreach my $date (sort keys %{$indHash->{$sp}->{$station}}) {
					foreach my $ind (sort keys %{$indHash->{$sp}->{$station}->{$date}}) {
						my $localHash =  \%{$indHash->{$sp}->{$station}->{$date}->{$ind}};
						if($localHash->{nmaster} == $nRequiredMaster && $localHash->{noptional} >= $minOptionalGenes ) {
							#this one passes the tests
							#calculate the overall ACGT score
							my $totalACGT = 0;
							foreach my $gene (@allGenes) {
								#####
	# 							print "$gene\n";
								if(exists $localHash->{$gene}) {
									my $seqid = $localHash->{$gene}->{id};
									$totalACGT += $localHash->{$gene}->{seqscore};
								}
							}
							if($totalACGT > $bestACGTScore) {
								$bestIndHash = \%{$indHash->{$sp}->{$station}->{$date}->{$ind}}
							}
						}
					}
				}

				#report the best individual if any
				if(defined $bestIndHash) {
					print TAB "$unitId";
					++$nIncluded;
					foreach my $gene (@allGenes) {
						if(exists $bestIndHash->{$gene}) {
							my $seqid = $bestIndHash->{$gene}->{id};
							print TAB "\t$seqid";
							$seqToExport{$gene}{$seqid} = 1;
						}
						else {
							print TAB "\t";
						}
					}
					print TAB "\n";
				}
				else {
					++$nRejected;
				}
			}

			elsif($criteria eq 'bestchimera') {
				#get the best sequences foreach gene
				#then checked that we pass the cutoffs and print
				my %selectedSeq;
				my $nmaster = 0;
				my $noptional = 0;
				foreach my $gene (@allGenes) {
					my $bestSeqId;
					my $bestSeqScore = 0;

					foreach my $date (sort keys %{$indHash->{$sp}->{$station}}) {
						foreach my $ind (sort keys %{$indHash->{$sp}->{$station}->{$date}}) {
							my $localHash =  \%{$indHash->{$sp}->{$station}->{$date}->{$ind}};
							if(exists $localHash->{$gene}) {
									my $seqid = $localHash->{$gene}->{id};
									my $score = $localHash->{$gene}->{seqscore};
									if($score > $bestSeqScore) {
										$bestSeqScore = $score;
										$bestSeqId = $seqid;
									}
							}
						}
					}

					if($bestSeqId) {
						$selectedSeq{$gene} = $bestSeqId;
						++$nmaster if($geneCategories{$gene} eq 'master');
						++$noptional if($geneCategories{$gene} eq 'optional');
					}
				}
				#do we have all the required genes?
				if($nmaster == $nRequiredMaster && $noptional >= $minOptionalGenes ) {
					print TAB "$sp|$station";
					++$nIncluded;
					foreach my $gene (@allGenes) {
						if(exists $selectedSeq{$gene}) {
							my $seqid = $selectedSeq{$gene};
							print TAB "\t$seqid";
							$seqToExport{$gene}{$seqid} = 1;
						}
						else {
							print TAB "\t";
						}
					}
					print TAB "\n";
				}
				else {
					++$nRejected;
				}
			}

		}
	}
}


elsif($unit eq 'species') {
	foreach my $sp (sort keys %$indHash) {
		#this is the unit level
# 		get all the individual below this level, and either select the best indiv or the best chimera
		my $unitId = $sp;
		#skip if a non selected unit
		if($selectedUnitsFile) {
			unless(map($unitId =~ /$_/, keys %unitToKeep)) { next }
		}
		if($criteria eq 'bestindiv') {
			my $bestIndHash;
			my $bestACGTScore = 0;
			foreach my $station (sort keys %{$indHash->{$sp}}) {
				foreach my $date (sort keys %{$indHash->{$sp}->{$station}}) {
					foreach my $ind (sort keys %{$indHash->{$sp}->{$station}->{$date}}) {
						my $localHash =  \%{$indHash->{$sp}->{$station}->{$date}->{$ind}};
						if($localHash->{nmaster} == $nRequiredMaster && $localHash->{noptional} >= $minOptionalGenes ) {
							#this one passes the tests
							#calculate the overall ACGT score
							my $totalACGT = 0;
							foreach my $gene (@allGenes) {
								#####
	# 							print "$gene\n";
								if(exists $localHash->{$gene}) {
									my $seqid = $localHash->{$gene}->{id};
									$totalACGT += $localHash->{$gene}->{seqscore};
								}
							}
							if($totalACGT > $bestACGTScore) {
								$bestIndHash = \%{$indHash->{$sp}->{$station}->{$date}->{$ind}}
							}
						}
					}
				}
			}

			#report the best individual/species if any
			if(defined $bestIndHash) {
				print TAB "$sp";
				++$nIncluded;
				foreach my $gene (@allGenes) {
					if(exists $bestIndHash->{$gene}) {
						my $seqid = $bestIndHash->{$gene}->{id};
						print TAB "\t$seqid";
						$seqToExport{$gene}{$seqid} = 1;
					}
					else {
						print TAB "\t";
					}
				}
				print TAB "\n";
			}
			else {
				++$nRejected;
			}
		}

		elsif($criteria eq 'bestchimera') {
			#get the best sequences foreach gene
			#then checked that we pass the cutoffs and print
			my %selectedSeq;
			my $nmaster = 0;
			my $noptional = 0;
			foreach my $gene (@allGenes) {
				my $bestSeqId;
				my $bestSeqScore = 0;

				foreach my $station (sort keys %{$indHash->{$sp}}) {
					foreach my $date (sort keys %{$indHash->{$sp}->{$station}}) {
						foreach my $ind (sort keys %{$indHash->{$sp}->{$station}->{$date}}) {
							my $localHash =  \%{$indHash->{$sp}->{$station}->{$date}->{$ind}};
							if(exists $localHash->{$gene}) {
									my $seqid = $localHash->{$gene}->{id};
									my $score = $localHash->{$gene}->{seqscore};
									if($score > $bestSeqScore) {
										$bestSeqScore = $score;
										$bestSeqId = $seqid;
									}
							}
						}
					}
				}

				if($bestSeqId) {
					$selectedSeq{$gene} = $bestSeqId;
					++$nmaster if($geneCategories{$gene} eq 'master');
					++$noptional if($geneCategories{$gene} eq 'optional');
				}
			}
			#do we have all the required genes?
			if($nmaster == $nRequiredMaster && $noptional >= $minOptionalGenes ) {
				print TAB "$sp";
				++$nIncluded;
				foreach my $gene (@allGenes) {
					if(exists $selectedSeq{$gene}) {
						my $seqid = $selectedSeq{$gene};
						print TAB "\t$seqid";
						$seqToExport{$gene}{$seqid} = 1;
					}
					else {
						print TAB "\t";
					}
				}
				print TAB "\n";
			}
			else {
				++$nRejected;
			}
		}
	}
}


elsif($unit eq 'motu') {
	#first get all the individuals by Motu and populate a new hash where species are replaced by motus
	my $motuHash;
	my @allmotus = sort keys %motu2seq;
	foreach my $sp (sort keys %$indHash) {
		foreach my $station (sort keys %{$indHash->{$sp}}) {
			foreach my $date (sort keys %{$indHash->{$sp}->{$station}}) {
				foreach my $ind (sort keys %{$indHash->{$sp}->{$station}->{$date}}) {
					my $localHash =  \%{$indHash->{$sp}->{$station}->{$date}->{$ind}};
					if(defined $localHash->{MOTU}) {
						my $motu = $localHash->{MOTU};
# 						push @{$motuHash{$motu}}, $localHash;
						$motuHash->{$motu}->{$station}->{$date}->{$ind} = $localHash;
					}
				}
			}
		}
	}
	open OUT4, ">motuHash.yaml";
	print OUT4 Dump($motuHash);
	#then work from this hash now on
	foreach my $motu (sort keys %$motuHash) {
		#this is the unit level
# 		get all the individual below this level, and either select the best indiv or the best chimera
		my $unitId = $motu;
# 		next if($selectedUnitsFile && !exists($unitToKeep{$unitId})); #does not make sense here, as we consider that the MOTU file is also used for this
		if($criteria eq 'bestindiv') {
			my $bestIndHash;
			my $bestACGTScore = 0;
			foreach my $station (sort keys %{$motuHash->{$motu}}) {
				foreach my $date (sort keys %{$motuHash->{$motu}->{$station}}) {
					foreach my $ind (sort keys %{$motuHash->{$motu}->{$station}->{$date}}) {
						my $localHash =  \%{$motuHash->{$motu}->{$station}->{$date}->{$ind}};
						if($localHash->{nmaster} == $nRequiredMaster && $localHash->{noptional} >= $minOptionalGenes ) {
							#this one passes the tests
							#calculate the overall ACGT score
							my $totalACGT = 0;
							foreach my $gene (@allGenes) {
								#####
	# 							print "$gene\n";
								if(exists $localHash->{$gene}) {
									my $seqid = $localHash->{$gene}->{id};
									$totalACGT += $localHash->{$gene}->{seqscore};
								}
							}
							if($totalACGT > $bestACGTScore) {
								$bestIndHash = \%{$motuHash->{$motu}->{$station}->{$date}->{$ind}}
							}
						}
					}
				}
			}
			#report the best individual/species if any
			if(defined $bestIndHash) {
				print TAB "$motu";
				++$nIncluded;
				foreach my $gene (@allGenes) {
					if(exists $bestIndHash->{$gene}) {
						my $seqid = $bestIndHash->{$gene}->{id};
						print TAB "\t$seqid";
						$seqToExport{$gene}{$seqid} = 1;
					}
					else {
						print TAB "\t";
					}
				}
				print TAB "\n";
			}
			else {
				++$nRejected;
			}
		}
		elsif($criteria eq 'bestchimera') {
			#get the best sequences foreach gene
			#then checked that we pass the cutoffs and print
			my %selectedSeq;
			my $nmaster = 0;
			my $noptional = 0;
			foreach my $gene (@allGenes) {
				my $bestSeqId;
				my $bestSeqScore = 0;
				foreach my $station (sort keys %{$motuHash->{$motu}}) {
					foreach my $date (sort keys %{$motuHash->{$motu}->{$station}}) {
						foreach my $ind (sort keys %{$motuHash->{$motu}->{$station}->{$date}}) {
							my $localHash =  \%{$motuHash->{$motu}->{$station}->{$date}->{$ind}};
							if(exists $localHash->{$gene}) {
								my $seqid = $localHash->{$gene}->{id};
								my $score = $localHash->{$gene}->{seqscore};
								if($score > $bestSeqScore) {
									$bestSeqScore = $score;
									$bestSeqId = $seqid;
								}
							}
						}
					}
				}
				if($bestSeqId) {
					$selectedSeq{$gene} = $bestSeqId;
					++$nmaster if($geneCategories{$gene} eq 'master');
					++$noptional if($geneCategories{$gene} eq 'optional');
				}
			}
			#do we have all the required genes?
			if($nmaster == $nRequiredMaster && $noptional >= $minOptionalGenes ) {
				print TAB "$motu";
				++$nIncluded;
				foreach my $gene (@allGenes) {
					if(exists $selectedSeq{$gene}) {
						my $seqid = $selectedSeq{$gene};
						print TAB "\t$seqid";
						$seqToExport{$gene}{$seqid} = 1;
					}
					else {
						print TAB "\t";
					}
				}
				print TAB "\n";
			}
			else {
				++$nRejected;
			}
		}
	}
}


#####exporting the sequences
if($exportseq) {
	foreach my $gene (keys %seqToExport) {
		my $geneFile = $geneToGeneFile{$gene};
		my $name = fileparse($geneFile);
		my $exportFile = $prefix . $name;
		my $seqin = Bio::SeqIO->new( -file => $geneFile, -format => 'fasta', -alphabet => 'dna');
		my $seqout = Bio::SeqIO->new( -file => ">$exportFile", -format => 'fasta');
		while(my $seq = $seqin->next_seq ) {
			my $id = $seq->id;
			$seqout->write_seq($seq) if(exists $seqToExport{$gene}{$id});
		}
	}
}

print "Number of units keept: $nIncluded
Number of units rejected: $nRejected\n";

#####SUB

sub calculate_seqscore {
	my $str = shift;
	my $n = 0;
	while ($str =~ /[ACTG]/gi) { ++$n };
	return $n;
}

sub calculate_ambiguityscore {
	my $str = shift;
	my $n = 0;
	while ($str =~ /[nkmbvswdyrh\?]/gi) { ++$n };
	return $n;
}



sub parse_formatv3 {
    my ($name) = @_;
    if($name =~ /^(.+)\|(\w+)_(\d+)_(\w+)_(\w+)\|(.*)$/) {
		my $head = $1;
		my $loc = $2;
		my $date = $3;
		my $ind = $4;
		my $first_yas = $5;
		my $tail = $6;
		my @yas;
		@yas = search_yas_in_tail($tail) unless (!defined $tail or $ind eq 'RNAseq');
		push @yas, $first_yas;

		#if this is a NCBI ref assign the AN to the ind code (we should fix this one day)
		if($ind eq 'NCBI' or $ind eq 'RNAseq') {
			$ind = $first_yas;
		}

		my $hash;
		$hash = {
			head => $head,
			loc => $loc,
			date => $date,
			ind => $ind,
			tail => $tail,
			original => $name,
			format => 'v3',
			firstyas => $first_yas,
			yas => \@yas,
		};
		return $hash;
	}
	else { print "This sequence does not fit the E3S format: $name\n"; }
}


sub parse_formatv4 {
    my ($name) = @_;
		    if($name =~ /^(.+)\|(\w+)_(\d+)_(\w+)_(\w+)\|(.*)$/) {
				my $head = $1;
				my $loc = $2;
				my $date = $3;
				my $ind = $4;
				my $first_yas = $5;
				my $tail = $6;
				my @yas;

				#search for sequences not produced locally by sanger seq
				if(!defined $tail or $tail =~ /^NCBI/i or $first_yas =~ /asb/) {
						my $seqType = 'NCBI';
						$seqType = 'transcript' if $first_yas =~ /asb/;
						if(!defined $tail) {
								$seqType = '?';
								print "Empty tail for $name\n" if !defined $tail;
						}

						if($seqType eq 'transcript') {
							$first_yas = $tail;
						}


						my $hash;
						$hash = {
							head => $head,
							loc => $loc,
							date => $date,
							ind => $ind,
							tail => $tail,
							original => $name,
							format => 'v4',
							firstyas => $first_yas,
							yas => \@yas,
							seqType => $seqType,
						};
						return $hash;

				} else {
							my $seqType = 'sanger';
							#decompose the tail
							@yas = search_yas_in_tail($tail);
							push @yas, $first_yas;
							#my $yas_primer = search_yas_primer_in_tail($first_yas, $tail, $name);

							my $hash;
							$hash = {
								head => $head,
								loc => $loc,
								date => $date,
								ind => $ind,
								tail => $tail,
								original => $name,
								format => 'v4',
								firstyas => $first_yas,
								yas => \@yas,
								#yas_primer => $yas_primer,
								seqType => $seqType,
							};
							return $hash;
			}
	}
	else {
		print "This sequence does not fit the E3S format: $name\n";
		return "error";
				}
}




sub search_yas_in_tail {
    my ($tail) = @_;
    my @yas;
    while($tail =~  /([YC][A-Z]{2}\d+)_/g) {
	push @yas, $1;
    }
    return(@yas);
}


__END__

=head1 NAME

chimera_assembler.pl

=head1 DESCRIPTION

Produce a sequence assembly table at the level of individuals,
populations, species or MOTUs

=head1 VERSION


v4 2017-june-16


=head1 SYNOPSIS

chimera_assembler.pl -master 16S_20120506.fas -master COI_20120506.fas -optional Opsine1_20120506.fas -unit station

=head1 OPTIONS

=over 8

=item B<-master>

<gene fasta file>, sequences that must be present

=item B<-optional>

<gene fasta file>, optional sequences

=item B<-min>

minimum number of sequences from the optional genes [0]

=item B<-unit>

<indiv|station|species|motu>, defines the unit level at which the assembly [indiv]

=item B<-criteria>

<bestindiv|bestchimera>, for any assembly unit higher than the
individual, either choose the the best individual or the best chimera
(ie. best set of sequences irrespective of the individual they belong to).
The decision is made based on ACGT counts [bestindiv]

=item B<-maxambiguity>

maximum percentage of ambiguity code a seq can have to be selected [5]

=item B<-useW>

use the sequences that have the working flag (W_). By default theses sequences are not used

=item B<-outtab>

 name of the output assembly [assembly_table.tsv]

=item B<-selectedunits>

name of a file listing the units to keep, the others will
discarded. The format should be the same as for the sequences (Species, Species|Station, ...)

=item B<-exportseq>

export fasta files with the sequences selected in the assembly table

=item B<-prefix>

prefix given to the output sequence files in addition to the original name [chimera_assembler_]

=item B<-motu>

<filename>, directs to a file that has at least two columns: 1) sequence name in
 E3S format (core or full), 2) MOTU name. If a single MOTU is found in a sp|location, the MOTU
propagated to the whole sp|location.

=item B<-noNCBI>

Do not use NCBI sequences

=item B<-noTranscript>

Do not use RNA-seq contigs

=back

=head1 AUTHOR

Tristan Lefebure, Tristan.Lefebure@univ-lyon1.fr

=cut
