#!/usr/bin/env perl
use strict;
use warnings;
######################################
#
# Authors : Aaron Darling and Guillaume Jospin
#
#
#
#
# demultiplex_dualBC.pl <illumina_directory> <barcode_list>
#
# Illumina_directory : Directory where the illumina files are located for the dual barcode demultiplex
#                      For a MiSeq run, files should look like :
#                      XTDM1_NoIndex_L001_R1_001.fastq.gz
#                      XTDM1_NoIndex_L001_R2_001.fastq.gz
#                      XTDM1_NoIndex_L001_R3_001.fastq.gz
#                      XTDM1_NoIndex_L001_R4_001.fastq.gz
#     FastQ files are assumed to be gzipped.  .gz filenames
#
#
# Barcode_list : list of the barcode name and barcode.  <barcode_label><tab><barcode>
#
# Printing to summary.txt the read counts for each barcode
# Mismatched_barcode : one of the two mates had a barcode that was not recognized from the list
# BC1.BC2 : Whenever a barcode pair was not listed in the name mapping file
#
#
######################################
## Use FLASH <read1> <read2> -m 30 -M 70 -x 0.25 -p 33 -o <merged file> -r 250 -f 250 -s 25
## Change -p 33 to whatever type calls for in the script.
use IO::Zlib;
use Getopt::Long;
my ( $full_print, $reverse, $fragment_length, $skip_merge );
my $read_length       = 300;
my $quality_threshold = 20;
my $fragment_std      = 20;
my $min_overlap       = 10;
my $max_overlap       = 70;
my $mismatch_ratio    = 0.25;
my $phred_value       = 33;
my $no_mismatch       = 0;
my $min_read_length   = 150;
my $skip_interleaved  = 0;
my $adaptor_file;
my $trim_path;
GetOptions(
			"full-print"        => \$full_print,
			"reverse"           => \$reverse,
			"q=i"               => \$quality_threshold,
			"read-len=i"        => \$read_length,
			"frag-len=i"        => \$fragment_length,
			"frag-std=i"        => \$fragment_std,
			"min-overlap=i"     => \$min_overlap,
			"max-overlap=i"     => \$max_overlap,
			"mismatch-ratio"    => \$mismatch_ratio,
			"skip-merge"        => \$skip_merge,
			"phred=i"           => \$phred_value,
			"min-read-length=i" => \$min_read_length,
			"no-mismatch"       => \$no_mismatch,
			"skip-interleaved"  => \$skip_interleaved,
    "trim-file=s" => \$adaptor_file,
    "trim-tool=s" => \$trim_path,
);
my $usage = "Wrong number of arguments\nUsage:\ndemultiplex_dualBC.pl <options> <illumina_directory> <mapping_file> <output_directory> <filename_core>\n";
die("$usage")                            if @ARGV != 4;
print STDERR "Not allowing mismatches\n" if $no_mismatch;
print STDERR "Indices are switched\n"    if $reverse;

die "$ARGV[1] was not found\n" unless -e $ARGV[1];

#reading barcodes
#my %barcode_forward = ();
my %barcode_rev   = ();
my $output_dir    = $ARGV[2];
my $out_file_core = $ARGV[3];
##if the $output_dir does not exists create it.
print STDERR "Creating $output_dir\n" unless -e $output_dir;
`mkdir -p $output_dir`                unless -e $output_dir;
##if the $output_dir does not exists create it.
print STDERR "Creating $output_dir/interleaved_fastq\n" unless -e "$output_dir/interleaved_fastq";
`mkdir -p $output_dir/interleaved_fastq`                unless -e "$output_dir/interleaved_fastq";
##if the $output_dir does not exists create it.
print STDERR "Creating $output_dir/raw_fastq\n" unless -e $output_dir."/raw_fastq";
`mkdir -p $output_dir/raw_fastq`                unless -e $output_dir."/raw_fastq";
##if the $output_dir does not exists create it.
print STDERR "Creating $output_dir/qiime_ready\n" unless -e $output_dir."/qiime_ready";
`mkdir -p $output_dir/qiime_ready`                unless -e $output_dir."/qiime_ready";

##if the $output_dir/trimmed_fastq does not exists create it.
print STDERR "Creating $output_dir/trimmed_fastq\n" unless -e $output_dir."/trimmed_fastq";
`mkdir -p $output_dir/trimmed_fastq`                unless -e $output_dir."/trimmed_fastq";

print STDERR "Reading barcodes and sample mapping\n";

open( INBC, $ARGV[1] );
my $header = <INBC>;
$header =~ s/^#//;
my @cols                    = split( /\s+/, $header );
my %mapping                 = ();
my %output_filehandles_1    = ();
my %output_filehandles_2    = ();
my %output_filehandles_full = ();
my %sample_read_count       = ();
my %mapping_file            = ();
my %total_codes             = ();
my %barcode                 = ();
my %samples                 = ();
my $line_count              = 0;
while (<INBC>) {
	chomp($_);
	next if $_ eq "";
	$line_count++;
	my @line = split( /\t/, $_ );
	for ( my $i = 0; $i < scalar(@line); $i++ ) {
		my $key   = $cols[$i];
		my $value = $line[$i];
		$mapping_file{ $line[0] }{$key} = $value;
	}
	unless (exists $mapping_file{ $line[0] }{"BarcodeSequence"} && exists $mapping_file{ $line[0] }{"ReverseBarcode"}){
	    print STDERR "**** WARNING **** Did not find BarcodeSequence or ReverseBarcode values.\nSkipping $_\n";
	    next;
	}
	$samples{ $line[0] } = 1;
	#changes spaces to . in sample IDs
	$line[0] =~ s/\s/\./g;
	if ( length $mapping_file{ $line[0] }{"BarcodeSequence"} != 8 ) {
		print STDERR "*** Warning *** Barcode sequence $mapping_file{ $line[0] }{\"BarcodeSequence\"} is not 8 bases long\n";
	}
	if ( length $mapping_file{ $line[0] }{"ReverseBarcode"} != 8 ) {
		print STDERR "*** Warning *** Barcode sequence $mapping_file{ $line[0] }{\"ReverseBarcode\"} is not 8 bases long\n";
	}
	if ( $mapping_file{ $line[0] }{"BarcodeSequence"} =~ m/[^ATCG]/ ) {
		die "*** Exiting *** Barcode sequence $mapping_file{ $line[0] }{\"BarcodeSequence\"} has illigal characters (not ATCG)\n";
	}
	if ( $mapping_file{ $line[0] }{"ReverseBarcode"} =~ m/[^ATCG]/ ) {
		die "*** Exiting *** Barcode sequence".$mapping_file{ $line[0] }{"ReverseBarcode"}." has illigal characters (not ATCG)\n";
	}
	my $code_original = $mapping_file{ $line[0] }{"BarcodeSequence"}.$mapping_file{ $line[0] }{"ReverseBarcode"};
	print STDERR "WARNING: $line[0]\t$code_original was already defined for $barcode{$code_original}\n"
	  if defined $barcode{$code_original} && $barcode{$code_original} ne $line[0];
	$barcode{$code_original} = $line[0];
	unless ($no_mismatch) {

		# insert all single-error barcodes
		for ( my $i = 0; $i < length( $mapping_file{ $line[0] }{"BarcodeSequence"} ); $i++ ) {
			my @chars = ( "A", "C", "G", "T", "N" );
			my $s = $mapping_file{ $line[0] }{"BarcodeSequence"};
			foreach my $ck (@chars) {
				substr( $s, $_, 1 ) =~ s/[ACGT]/$ck/ for $i;
				for ( my $j = 0; $j < length( $mapping_file{ $line[0] }{"ReverseBarcode"} ); $j++ ) {
					my @chars2 = ( "A", "C", "G", "T", "N" );
					my $s2 = $mapping_file{ $line[0] }{"ReverseBarcode"};
					foreach my $ck2 (@chars2) {
						substr( $s2, $_, 1 ) =~ s/[ACGT]/$ck2/ for $j;
						my $long_code = $s.$s2;
						print STDERR "WARNING: $line[0]\t$long_code was already defined for $barcode{$long_code}\n"
						  if defined $barcode{$long_code} && $barcode{$long_code} ne $line[0];
						$barcode{$long_code} = $line[0];
					}
				}
			}
		}
	}
}
print STDERR "Found ".scalar( keys %samples )." samples. Read $line_count non-empty lines in the mapping file\n";
print STDERR "*** Warning ***  Lots of samples within the same mapping file may not work.  A safe batch of samples is 96 or lower.\n" if $line_count > 96;
print STDERR "Total number of barcodes : ".scalar( keys %barcode );
print STDERR " using 1 maximum mismatch per barcode" unless $no_mismatch;
print STDERR "\n";
close(INBC);

#foreach my $fw ( keys %mapping ) {
#	foreach my $rv ( keys %{ $mapping{$fw} } ) {
#		print "$fw\t$rv\t$mapping{$fw}{$rv}\n";
#	}
#}
my %bc_count = ();
# Opening filehandles for future print statements
foreach my $code ( sort keys %barcode ) {
	unless ($skip_interleaved) {
		$output_filehandles_full{ $barcode{$code} } = new IO::Zlib;
		$output_filehandles_full{ $barcode{$code} }->open( "$output_dir/interleaved_fastq/$out_file_core"."_$barcode{$code}.fastq.gz", "wb9" );
	}
	open( $output_filehandles_1{ $barcode{$code} }, ">$output_dir/raw_fastq/$out_file_core"."_$barcode{$code}"."_1".".fastq" ) unless $skip_merge;
	open( $output_filehandles_2{ $barcode{$code} }, ">$output_dir/raw_fastq/$out_file_core"."_$barcode{$code}"."_2".".fastq" ) unless $skip_merge;
	$bc_count{$barcode{$code}}=0;
}

## FastQ files are assumed to be gzipped.  .gz filenames

my @files    = <$ARGV[0]/*_R1_*.fastq.gz>;
foreach my $file (@files) {
	print STDERR "Processing $file\n";
	$file =~ m/^(\S+)_\S\S_(\d+).fastq.gz/;
	my $core  = $1;
	my $index = $2;
	open( my $TYPETEST, "zcat $file |" );
	my $type = get_sequence_input_type($TYPETEST);

	#print STDERR "TYPE :".$type->{qtype}."\n";
	close($TYPETEST);
	my $read_type = $type->{qtype};
	open( READ1,  "zcat $file |" );
	open( READ2,  "zcat $core"."_R4_$index.fastq.gz |" );
	open( INDEX1, "zcat $core"."_R2_$index.fastq.gz |" );
	open( INDEX2, "zcat $core"."_R3_$index.fastq.gz |" );
	my @read1  = ();
	my @read2  = ();
	my @index1 = ();
	my @index2 = ();
	## read the reads and indices from the input files.
	while (1) {
		for ( my $i = 0; $i < 4; $i++ ) {
			if ($reverse) {
				$read1[$i]  = <READ2>;
				$read2[$i]  = <READ1>;
				$index1[$i] = <INDEX2>;
				$index2[$i] = <INDEX1>;
			} else {
				$read1[$i]  = <READ1>;
				$read2[$i]  = <READ2>;
				$index1[$i] = <INDEX1>;
				$index2[$i] = <INDEX2>;
			}
		}
		## stop reading if we reached the end of the files
		last if !defined( $read1[0] );
		## figure out what barcodes we are dealing with
		my $i1 = $index1[1];
		my $i2 = $index2[1];
		chomp($i1);
		chomp($i2);
		$total_codes{$i1}{$i2} = 0 unless exists $total_codes{$i1}{$i2};
		$total_codes{$i1}{$i2}++;
		my $long_code = $i1.$i2;

		if ( exists $barcode{$long_code} ) {

			#print "FOUND CODE\n";
			$bc_count{ $barcode{$long_code} } = 0 unless exists $bc_count{ $barcode{$long_code} };
			$bc_count{ $barcode{$long_code} }++;
			my $READ_HANDLE_1    = $output_filehandles_1{ $barcode{$long_code} };
			my $READ_HANDLE_2    = $output_filehandles_2{ $barcode{$long_code} };
			my $READ_HANDLE_FULL = $output_filehandles_full{ $barcode{$long_code} } unless $skip_interleaved;
			## quality trim the reads
			qtrim_read( read => \@read1, quality => $quality_threshold, readtype => $type );
			qtrim_read( read => \@read2, quality => $quality_threshold, readtype => $type );
			## print the trimmed reads to an intermediate file if specified in $full_print

			## merge the reads if possible
			## my @merged_read = align_and_merge_reads(read1=> \@read1, read2=> \@read2 );

			## print the merged reads to an intermediate file if specified in $full_print
			## discard the reads if it wasn't possible to merge
			$read1[0] = clean_line( line => $read1[0], num => 1 );
			$read2[0] = clean_line( line => $read2[0], num => 2 );

			## Change the read header to accomodate for barcoding.
			## print the reads to their respective files
			print $READ_HANDLE_FULL @read1 if exists $output_filehandles_full{ $barcode{$long_code} } && !$skip_interleaved;
			print $READ_HANDLE_FULL @read2 if exists $output_filehandles_full{ $barcode{$long_code} } && !$skip_interleaved;
			print $READ_HANDLE_1 @read1    if exists $output_filehandles_1{ $barcode{$long_code} };
			print $READ_HANDLE_2 @read2    if exists $output_filehandles_2{ $barcode{$long_code} };
		} else {
			$bc_count{mismatch} = 0 unless exists $bc_count{mismatch};
			$bc_count{mismatch}++;
		}
	}
}
##close files and flush IO buffers
foreach my $handle ( keys(%output_filehandles_1) ) {
	$output_filehandles_1{$handle}->close();
	$output_filehandles_2{$handle}->close();
	if ( -z "$output_dir/raw_fastq/$out_file_core"."_$handle"."_1.fastq" ) {
		print STDERR "$output_dir/raw_fastq/$out_file_core"."_$handle"."_1.fastq"." is empty, can't align and merge reads\n";
		next;
	}
	if ( -z "$output_dir/raw_fastq/$out_file_core"."_$handle"."_2.fastq" ) {
		print STDERR "$output_dir/raw_fastq/$out_file_core"."_$handle"."_2.fastq"." is empty, can't align and merge reads\n";
		next;
	}

	my $source_dir= "raw_fastq";
	if($trim_path && $adaptor_file){
	    # Need to trim out adaptors in case we read through the other end of the read.
	    # Trimming the content of $adaptor_file
	    my $trim_cmd = "java -jar $trim_path/trimmomatic.jar PE -phred33 -baseout \"$output_dir/trimmed_fastq/$out_file_core"."_$handle\" \"$output_dir/raw_fastq/$out_file_core"."_$handle"."_1.fastq\" \"$output_dir/raw_fastq/$out_file_core"."_$handle"."_2.fastq\" ILLUMINACLIP:$adaptor_file:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 >/dev/null 2>/dev/null";
	    `$trim_cmd`;
	    $source_dir = "trimmed_fastq";
	}

	
	
	#print "HANDLE : $handle\n";
	#print "PHRED VALUE : $phred_value\n";
	my $options   = "-m $min_overlap -M $max_overlap -p $phred_value -s $fragment_std -r $read_length -x $mismatch_ratio";
	my $flash_cmd ="flash \"$output_dir/$source_dir/$out_file_core"."_$handle";
	$flash_cmd .= $source_dir eq "trimmed_fastq" ? "_1P\" ":"_1.fastq\" ";
#	  ."_1P" if ($source_dir eq "trimmed_fastq");
#	  ."_1.fastq\"";
	$flash_cmd .= "\"$output_dir/$source_dir/$out_file_core"."_$handle";
	$flash_cmd .= $source_dir eq "trimmed_fastq" ? "_2P\" ":"_2.fastq\" ";
#	  ."_2P" if ($source_dir eq "trimmed_fastq");
#	  ."_2.fastq\""; 
	$flash_cmd .="$options -d \"$output_dir/qiime_ready\" -o \"$out_file_core"
	  .".$handle\" > /dev/null 2> /dev/null";

	print STDERR "RUNNING : $flash_cmd\n";
	system($flash_cmd)
	  unless -z "$output_dir/raw_fastq/$out_file_core"."_$handle"."_1.fastq" && -z "$output_dir/raw_fastq/$out_file_core"."_$handle"."_2.fastq";
	open( OUT_QIIME,     ">$output_dir/qiime_ready/$out_file_core.M.".$handle.".faa" );
	open( OUT_QIIME_FOR, ">$output_dir/qiime_ready/$out_file_core.MPF.".$handle.".faa" );
	open( OUT_QIIME_REV, ">$output_dir/qiime_ready/$out_file_core.MPR.".$handle.".faa" );
	
	$sample_read_count{$handle} = 0;
	if ( -e "$output_dir/qiime_ready/$out_file_core.".$handle.".extendedFrags.fastq" ) {
		open( INMERGED, "$output_dir/qiime_ready/$out_file_core.".$handle.".extendedFrags.fastq" );
		my @read = ();
		while (<INMERGED>) {
			$read[0] = $_;
			$read[1] = <INMERGED>;
			$read[2] = <INMERGED>;
			$read[3] = <INMERGED>;
			my @qiime_read1 = convert_to_qiime_read( read_array => \@read, sample => $handle );
			$qiime_read1[0] =~ s/\n/\tMerged\n/;

			#		print STDERR "qiime_Ready ready @qiime_read1";
			print OUT_QIIME_REV @qiime_read1;
			print OUT_QIIME_FOR @qiime_read1;
			print OUT_QIIME @qiime_read1;
			@read = ();
		}
		close(INMERGED);
	} else {
		print STDERR "$output_dir/qiime_ready/$out_file_core.".$handle.".extendedFrags.fastq was not found\n";
		next;
	}
	my $count = 0;
	if ( -e "$output_dir/qiime_ready/$out_file_core.".$handle.".notCombined_1.fastq" ) {
		open( INFORWARD, "$output_dir/qiime_ready/$out_file_core.".$handle.".notCombined_1.fastq" );
		my @read = ();
		while (<INFORWARD>) {
			$read[0] = $_;
			$read[1] = <INFORWARD>;
			$read[2] = <INFORWARD>;
			$read[3] = <INFORWARD>;

			#		print "Read $read[1]";
			qtrim_read( read => \@read, readtype => "phred".$phred_value, quality => $quality_threshold );

			#		print "Read $read[1]"."Lengh of new_Read = ".length( $read[1] )."\n\n\n";
			next if length( $read[1] ) < $min_read_length;
			my @qiime_read1 = convert_to_qiime_read( read_array => \@read, sample => $handle );

			#print STDERR "qiime_Ready ready @qiime_read1";
			$qiime_read1[0] =~ s/\n/\tNot_merged\n/;
			print OUT_QIIME_FOR @qiime_read1;
			@read = ();
			$count++;
		}
		close(INFORWARD);
	} else {
		print STDERR "$output_dir/qiime_ready/$out_file_core.".$handle.".notCombined_1.fastq was not found\n";
	}
	##reset the read count increment
	$sample_read_count{$handle} = $sample_read_count{$handle} - $count;

	if ( -e "$output_dir/qiime_ready/$out_file_core.".$handle.".notCombined_2.fastq" ) {
		open( INREVERSE, "$output_dir/qiime_ready/$out_file_core.".$handle.".notCombined_2.fastq" );
		my @read = ();
		while (<INREVERSE>) {
			$read[0] = $_;
			$read[1] = <INREVERSE>;
			$read[2] = <INREVERSE>;
			$read[3] = <INREVERSE>;
			qtrim_read( read => \@read, readtype => "phred".$phred_value, quality => $quality_threshold );
			next if length( $read[1] ) < $min_read_length;
			my @qiime_read1 = convert_to_qiime_read( read_array => \@read, sample => $handle );

			#print STDERR "qiime_Ready ready @qiime_read1";
			$qiime_read1[0] =~ s/\n/\tNot_merged\n/;
			print OUT_QIIME_REV @qiime_read1;
			@read = ();
			$count++;
		}
		close(INREVERSE);
	} else {
		print STDERR "$output_dir/qiime_ready/$out_file_core.".$handle.".notCombined_1.fastq was not found\n";
	}
	close(OUT_QIIME);

	#zip the raw files
	my $file_1 = "$output_dir/raw_fastq/$out_file_core"."_$handle"."_1.fastq"; 
	my $file_2 = "$output_dir/raw_fastq/$out_file_core"."_$handle"."_2.fastq"; 
	`gzip -f "$file_1"`;
	`gzip -f "$file_2"`;
}

## Write a summary of the run
open( OUT, ">$output_dir/summary.txt" );
foreach my $code ( sort { $bc_count{$a} <=>$bc_count{$b} } keys %bc_count ) {
	print OUT $code."\t".$bc_count{$code}."\n";
}
close(OUT);

open( RUNINFO, ">$output_dir/run_info.txt" );
print RUNINFO "Read length : $read_length\n";
print RUNINFO "Read quality threshold : $quality_threshold\n";
print RUNINFO "fragment standard deviation : $fragment_std\n";
print RUNINFO "Minimum overlap : $min_overlap\n";
print RUNINFO "Maximum overlap : $max_overlap\n";
print RUNINFO "Mismatch ratio : $mismatch_ratio\n";
print RUNINFO "Phred score : PHRED$phred_value\n";
print RUNINFO "Mismatches allowed : $no_mismatch\n";
print RUNINFO "Minimum read length : $min_read_length\n";
$reverse = 0 unless $reverse;
print RUNINFO "Reversed read : $reverse\n";
$skip_merge = 0 unless $skip_merge;
print RUNINFO "Skip merge : $skip_merge\n";
print RUNINFO "\n\n";
print RUNINFO "Read $line_count non-empty lines from the mapping file\n";
close(RUNINFO);

sub convert_to_qiime_read {
	my %args   = @_;
	my @read   = @{ $args{read_array} };
	my $sample = $args{sample};
	my @return_array;
	my $read_seq = $read[1];
	chomp($read_seq);
	##$read_seq = reverse($read_seq);
	##$read_seq =~ tr/ACGTacgt/TGCAtgca/;
	$sample_read_count{$sample} = 0 unless exists $sample_read_count{$sample};
	$sample_read_count{$sample}++;
	my $old_sample = $sample;
	$sample =~ s/_/./g;    #qiime does not like _ in the sample names.
	$return_array[0] = ">$sample"."_$sample_read_count{$old_sample}\n";
	$return_array[1] = $read_seq."\n";

	#print "@read\n@return_array";

	return @return_array;
}

sub clean_line {
	my %args = @_;
	$args{line} =~ s/ /:/g;
	chomp $args{line};
	$args{line} .= "/".$args{num}."\n";
	my @line = split( /:/, $args{line} );
	splice( @line, 7, 1 );
	$args{line} = join( ':', @line );
	return $args{line};
}

=head2 qtrim_read

trims a fastq read to a particular quality score using Heng Li's algorithm from bwa.
code based on SGA's implementation.

=cut

sub qtrim_read {
	my %args     = @_;
	my $read     = $args{read};
	my $q        = $args{quality};
	my $readtype = $args{readtype};

	$q += 33 if $readtype eq "phred33";
	$q += 64 if $readtype eq "phred64";

	# Perform a soft-clipping of the sequence by removing low quality bases from the
	# 3' end using Heng Li's algorithm from bwa

	my $seq = @$read[1];
	chomp $seq;
	my $qq = @$read[3];
	chomp $qq;
	my @qual          = split( //, $qq );
	my $endpoint      = 0;                  # not inclusive
	my $max           = 0;
	my $i             = length($seq) - 1;
	my $terminalScore = ord( $qual[$i] );

	# Only perform soft-clipping if the last base has qual less than $q
	return if ( $terminalScore >= $q );

	my $subSum = 0;
	while ( $i >= 0 ) {
		my $ps    = ord( $qual[$i] );
		my $score = $q - $ps;
		$subSum += $score;
		if ( $subSum > $max ) {
			$max      = $subSum;
			$endpoint = $i;
		}
		$i--;
	}

	# Clip the read
	@$read[1] = substr( $seq, 0, $endpoint )."\n";
	@$read[3] = substr( @$read[3], 0, $endpoint )."\n";
}

=head2 get_sequence_input_type

Checks whether input is FastA, FastQ, which quality type (33 or 64), and DNA or AA
Returns a hash reference with the values 'seqtype', 'format', and 'qtype' populated.

=cut

sub get_sequence_input_type {
	my $FILE = shift;
	my %type;
	my $counter    = 0;
	my $maxfound   = 0;
	my $dnacount   = 0;
	my $line_count = 0;
	$type{seqtype} = "dna";
	$type{format}  = "unknown";
	$type{qtype}   = "none";
	my $allcount = 0;
	my $sequence = 1;
	my $minq     = 255;    # minimum fastq quality score (for detecting phred33/phred64)
	my @lines;

	while ( my $line = <$FILE> ) {

		#print "$line\n";
		if ( $line =~ /^>/ ) {
			$maxfound = $counter > $maxfound ? $counter : $maxfound;
			$counter = 0;
			$type{format} = "fasta" if $type{format} eq "unknown";
		} elsif ( $line =~ /^@/ || $line =~ /^\+/ ) {
			$counter  = 0;
			$sequence = 1;
			$type{format} = "fastq" if $type{format} eq "unknown";

			#print STDERR "File format detected is fastq\n";
		} elsif ( $line =~ /^\+/ ) {
			$sequence = 0;
			$type{format} = "fastq" if $type{format} eq "unknown";
		} elsif ( $type{format} eq "fastq" && !$sequence ) {

			#print STDERR "Checking quality scores\n";

			# check whether qualities are phred33 or phred64
			for my $q ( split( //, $line ) ) {
				$minq = ord($q) if ord($q) < $minq;
			}
		} elsif ($sequence) {

			#$sequence =~ s/[-\.]//g;    #removing gaps from the sequences
			$line =~ s/[-\.]//g;
			$counter  += length($line) - 1;
			$dnacount += $line =~ tr/[ACGTUNacgtun]//;
			$allcount += length($line) - 1;
		}
		$line_count++;
		last if ( $line_count > 10 );
		last if ( $counter > 100000 );
	}

	$maxfound = $counter > $maxfound ? $counter : $maxfound;
	$type{seqtype} = "protein" if ( $dnacount < $allcount * 0.75 );
	$type{seqtype} = "dna"
	  if ( $type{format} eq "fastq" );    # nobody using protein fastq (yet)
	$type{qtype} = "phred64" if $minq < 255;
	$type{qtype} = "phred33" if $minq < 64;
	$type{paired} = 0;                    # TODO: detect interleaved read pairing
	return \%type;
}
